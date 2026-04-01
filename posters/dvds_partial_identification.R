# ============================================================
#  Partial Identification of E[Y(1) | X] under Hidden Confounding
#  DVDS Method: Dorn, Guo & Kallus (2023)
#  "Doubly-Valid/Doubly-Sharp Sensitivity Analysis for Causal
#   Inference with Unmeasured Confounding"
#
#  Sensitivity model: Marginal Sensitivity Model (MSM)
#  with parameter Lambda (Λ).
#
#  Under MSM(Λ), the true propensity odds are within Λ of the
#  nominal (estimated) propensity odds:
#
#    1/Λ ≤ [π(x)/(1-π(x))] / [ê(x)/(1-ê(x))] ≤ Λ
#
#  where π(x) is the TRUE propensity (unknown), ê(x) is the
#  NOMINAL propensity estimated from observed data.
#  Λ = 1  →  no hidden confounding (standard IPW).
#  Λ > 1  →  allows for increasing amounts of confounding.
# ============================================================

library(quantreg)   # quantile regression for DVDS augmentation
library(ggplot2)
library(dplyr)

set.seed(2024)
n <- 2000

# ------------------------------------------------------------
# 1.  Data Generating Process (DGP)
# ------------------------------------------------------------
#
#  X  ~ N(0,1)          observed covariate
#  U  ~ N(0,1)          HIDDEN confounder (analyst does not see U)
#
#  True propensity: P(A=1 | X, U) = logistic(0.5*X + 0.9*U)
#  → U is a strong confounder of A → Y
#
#  Potential outcomes:
#    Y(1) = 2 + X + 0.8*U + ε1     (treated)
#    Y(0) = 0.5*X        + ε0     (control)
#
#  Observed: Y = A·Y(1) + (1-A)·Y(0)

X  <- rnorm(n)
U  <- rnorm(n)                                  # HIDDEN — never used in estimation

true_ps <- plogis(0.5 * X + 0.9 * U)           # depends on U
A  <- rbinom(n, 1, true_ps)

Y1_pot <- 2 + X + 0.8 * U + rnorm(n, sd = 0.5) # potential outcome Y(1)
Y0_pot <- 0.5 * X          + rnorm(n, sd = 0.5) # potential outcome Y(0)
Y  <- A * Y1_pot + (1 - A) * Y0_pot             # observed outcome

cat("============================================================\n")
cat("  ORACLE (uses U — unobservable in practice)\n")
cat("============================================================\n")
cat(sprintf("  True E[Y(1)]       = %.3f\n", mean(Y1_pot)))
cat(sprintf("  True E[Y(0)]       = %.3f\n", mean(Y0_pot)))
cat(sprintf("  True ATE           = %.3f\n\n", mean(Y1_pot - Y0_pot)))

# ------------------------------------------------------------
# 2.  Analyst's world: only (X, A, Y) are observed
#     Estimate nominal propensity score ignoring U
# ------------------------------------------------------------

ps_model <- glm(A ~ X, family = binomial)       # misspecified: omits U
pi_hat   <- fitted(ps_model)
pi_hat   <- pmax(pi_hat, 0.02)                  # trim for numerical stability

naive_ipw <- mean(A * Y / pi_hat)               # Horwitz-Thompson estimator

cat("============================================================\n")
cat("  NAIVE ANALYSIS (ignores hidden confounding)\n")
cat("============================================================\n")
cat(sprintf("  Naive IPW E[Y(1)]  = %.3f  (biased due to hidden U)\n\n",
            naive_ipw))

# ------------------------------------------------------------
# 3.  DVDS Marginal Bounds on E[Y(1)]
# ------------------------------------------------------------
#
#  Under MSM(Λ), the SHARP identified set for E[Y(1)] is
#  [L(Λ), U(Λ)] where the bounds come from a doubly-valid
#  augmented estimator (Dorn, Guo & Kallus 2023, Sec. 3):
#
#  Key idea: choose a quantile level τ = Λ/(1+Λ) (upper bound)
#  or τ = 1/(1+Λ) (lower bound) for the treated-unit outcome
#  distribution.  Then the DVDS estimator is:
#
#    ψ̂_upper = (1/n) Σ_i [ q̂_{τ}(X_i) + w_i^up · A_i · (Y_i − q̂_{τ}(X_i)) ]
#    ψ̂_lower = (1/n) Σ_i [ q̂_{τ'}(X_i) + w_i^lo · A_i · (Y_i − q̂_{τ'}(X_i)) ]
#
#  where the MSM tilted weights are:
#    w_i^up = Λ / (Λ · π̂(X_i) + 1 − π̂(X_i))
#    w_i^lo = 1 / (π̂(X_i) + Λ · (1 − π̂(X_i)))
#
#  and q̂_τ(x) is a quantile regression of Y on X at level τ,
#  fit on TREATED units only.
#
#  "Doubly valid" means the bound is valid if EITHER:
#   (a) the propensity score model is correct, OR
#   (b) the outcome quantile model is correct.
#
#  "Sharp" means no tighter bound is possible under MSM(Λ).

dvds_marginal_bounds <- function(Y, A, X, pi_hat, Lambda) {

  tau_up  <- Lambda / (1 + Lambda)   # quantile for upper bound
  tau_lo  <- 1      / (1 + Lambda)   # quantile for lower bound

  # MSM tilted IPW weights (one per observation)
  w_up <- Lambda / (Lambda * pi_hat + (1 - pi_hat))
  w_lo <- 1      / (pi_hat + Lambda * (1 - pi_hat))

  # Quantile regression on treated units only
  df_treated <- data.frame(Y = Y[A == 1], X = X[A == 1])

  qfit_up <- rq(Y ~ X, tau = tau_up, data = df_treated)
  qfit_lo <- rq(Y ~ X, tau = tau_lo, data = df_treated)

  q_up <- predict(qfit_up, newdata = data.frame(X = X))  # predict for ALL units
  q_lo <- predict(qfit_lo, newdata = data.frame(X = X))

  # DVDS augmented estimators (IPW + quantile control variate)
  psi_up <- mean(q_up + w_up * A * (Y - q_up))
  psi_lo <- mean(q_lo + w_lo * A * (Y - q_lo))

  return(list(lower   = psi_lo,
              upper   = psi_up,
              width   = psi_up - psi_lo,
              q_up    = q_up,
              q_lo    = q_lo,
              w_up    = w_up,
              w_lo    = w_lo,
              tau_up  = tau_up,
              tau_lo  = tau_lo))
}

# Compute bounds over a sequence of Lambda values
Lambda_seq  <- c(1.0, 1.5, 2.0, 2.5, 3.0)
true_EY1    <- mean(Y1_pot)

marginal_results <- do.call(rbind, lapply(Lambda_seq, function(L) {
  b <- dvds_marginal_bounds(Y, A, X, pi_hat, L)
  data.frame(Lambda      = L,
             Lower       = b$lower,
             Upper       = b$upper,
             Width       = b$width,
             Contains_truth = (b$lower <= true_EY1 & true_EY1 <= b$upper))
}))

cat("============================================================\n")
cat(sprintf("  DVDS MARGINAL BOUNDS on E[Y(1)]  (true = %.3f)\n", true_EY1))
cat("============================================================\n")
print(marginal_results, row.names = FALSE, digits = 3)
cat("\n")

# ------------------------------------------------------------
# 4.  Conditional Bounds: E[Y(1) | X = x]  for each Lambda
# ------------------------------------------------------------
#
#  We evaluate the DVDS bound locally at each grid point x0
#  using a Nadaraya-Watson kernel smoother with Gaussian kernel.
#  This gives the identified SET for the conditional mean
#  potential outcome as a function of X.

dvds_conditional_bounds <- function(Y, A, X, pi_hat, Lambda,
                                    x_grid, bandwidth = 0.4) {

  tau_up <- Lambda / (1 + Lambda)
  tau_lo <- 1      / (1 + Lambda)

  w_up <- Lambda / (Lambda * pi_hat + (1 - pi_hat))
  w_lo <- 1      / (pi_hat + Lambda * (1 - pi_hat))

  # Global quantile regression on treated (extrapolates to full X range)
  df_tr  <- data.frame(Y = Y[A == 1], X = X[A == 1])
  q_up   <- predict(rq(Y ~ X, tau = tau_up, data = df_tr),
                    newdata = data.frame(X = X))
  q_lo   <- predict(rq(Y ~ X, tau = tau_lo, data = df_tr),
                    newdata = data.frame(X = X))

  # Local DVDS bounds at each grid point using kernel weighting
  res <- vapply(x_grid, function(x0) {

    kw <- dnorm((X - x0) / bandwidth)   # Gaussian kernel weights
    kw_norm <- kw / sum(kw)             # normalise to sum to 1

    upper <- sum(kw_norm * (q_up + w_up * A * (Y - q_up)))
    lower <- sum(kw_norm * (q_lo + w_lo * A * (Y - q_lo)))

    # Oracle: true E[Y(1) | X ≈ x0]
    true_cond <- sum(kw_norm * Y1_pot)

    # Naive: sample mean of Y among treated near x0
    kw_tr <- kw * A
    naive  <- if (sum(kw_tr) > 1e-8) sum(kw_tr * Y) / sum(kw_tr)
              else NA_real_

    c(lower = lower, upper = upper,
      true  = true_cond, naive = naive)
  }, FUN.VALUE = numeric(4))

  data.frame(X      = x_grid,
             Lower  = res["lower", ],
             Upper  = res["upper", ],
             True   = res["true",  ],
             Naive  = res["naive", ],
             Lambda = Lambda)
}

x_grid <- seq(-2.5, 2.5, length.out = 60)

cond_list <- lapply(Lambda_seq, function(L) {
  dvds_conditional_bounds(Y, A, X, pi_hat, L, x_grid, bandwidth = 0.4)
})
cond_df <- bind_rows(cond_list)
cond_df$Lambda_label <- factor(paste0("Λ = ", cond_df$Lambda))

# Extract reference curves (same for all Lambda)
ref <- cond_df[cond_df$Lambda == 1, ]

# ------------------------------------------------------------
# 5.  Plot 1 — Conditional identified sets for each Lambda
# ------------------------------------------------------------

p1 <- ggplot(cond_df, aes(x = X)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = Lambda_label),
              alpha = 0.28) +
  geom_line(aes(y = (Lower + Upper) / 2, color = Lambda_label),
            linewidth = 0.65, linetype = "dotted") +
  geom_line(data = ref,
            aes(x = X, y = True),
            color = "black", linewidth = 1.1,
            linetype = "dashed", inherit.aes = FALSE) +
  geom_line(data = ref,
            aes(x = X, y = Naive),
            color = "steelblue", linewidth = 0.9,
            inherit.aes = FALSE) +
  facet_wrap(~Lambda_label, ncol = 3) +
  scale_fill_manual(
    values = c("#FEE08B","#FDAE61","#F46D43","#D73027","#A50026"),
    guide  = "none") +
  scale_color_manual(
    values = c("#FEE08B","#FDAE61","#F46D43","#D73027","#A50026"),
    guide  = "none") +
  labs(
    title    = "DVDS Partial Identification of E[Y(1) | X] under Hidden Confounding",
    subtitle = paste0("Shaded = identified set under MSM(Λ)  ·  ",
                      "Dashed black = true E[Y(1)|X] (oracle)  ·  ",
                      "Blue = naive mean of Y among treated"),
    x        = "Observed covariate  X",
    y        = "E[Y(1) | X]",
    caption  = "DVDS estimator — Dorn, Guo & Kallus (2023)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text    = element_text(face = "bold", size = 12),
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

# ------------------------------------------------------------
# 6.  Plot 2 — Marginal bound width as a function of Lambda
# ------------------------------------------------------------

p2 <- ggplot(marginal_results, aes(x = Lambda)) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "#FDAE61", alpha = 0.4) +
  geom_line(aes(y = Lower), color = "#D73027", linewidth = 1) +
  geom_line(aes(y = Upper), color = "#D73027", linewidth = 1) +
  geom_hline(yintercept = true_EY1,
             color = "black", linetype = "dashed", linewidth = 0.9) +
  geom_hline(yintercept = naive_ipw,
             color = "steelblue", linetype = "solid", linewidth = 0.9) +
  annotate("text", x = 1.05, y = true_EY1 + 0.06,
           label = "True E[Y(1)]", hjust = 0, size = 3.5, color = "black") +
  annotate("text", x = 1.05, y = naive_ipw - 0.1,
           label = "Naive IPW", hjust = 0, size = 3.5, color = "steelblue") +
  scale_x_continuous(breaks = Lambda_seq) +
  labs(
    title    = "Marginal Identified Set for E[Y(1)] vs Sensitivity Parameter Λ",
    subtitle = "Shaded = [lower, upper] DVDS bound · Λ=1 means no hidden confounding",
    x        = "Sensitivity parameter  Λ",
    y        = "Identified set for E[Y(1)]",
    caption  = "As Λ increases, the identified set widens, reflecting greater uncertainty."
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(color = "grey40", size = 10)
  )

# Print plots
print(p1)
print(p2)

cat("============================================================\n")
cat("  INTERPRETATION GUIDE\n")
cat("============================================================\n")
cat("  Λ = 1 : assumes no hidden confounding (point-identified).\n")
cat("  Λ > 1 : allows unmeasured confounders to shift treatment\n")
cat("          odds by up to factor Λ.\n\n")
cat("  The DVDS bound is:\n")
cat("  • Sharp  — no tighter interval is possible under MSM(Λ).\n")
cat("  • Doubly valid — consistent if EITHER the propensity\n")
cat("    model OR the outcome quantile model is correctly\n")
cat("    specified (but not necessarily both).\n\n")
cat("  In this simulation, the analyst's propensity model is\n")
cat("  MISSPECIFIED (omits U).  The naive IPW estimate is biased\n")
cat("  downward.  The true E[Y(1)] falls inside the identified\n")
cat("  set once Λ is large enough to cover the true confounding.\n")
