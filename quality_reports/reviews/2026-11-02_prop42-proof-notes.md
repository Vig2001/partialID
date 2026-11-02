# Notes on Brieuc's shorter proof of Proposition 4.2

Source: `prop42_proof.pdf` (B. Lehmann, 29 Oct 2026), the three-step replacement for the
eight-step proof of `prop:lower-asvalid` in `partialID_Overleaf/Methods/derivations.tex`.
Discussion from 2 Nov 2026. Nothing in the paper has been changed yet.

Notation as in the note: $Z_k = \sqrt{n_k}(\hat L_k - L_k)/\sigma_k$ with cdf $F_k$,
$\varepsilon_k = \sup_x |F_k(x) - \Phi(x)|$, $T_n(\omega) = aZ_1 + bZ_2$ with $a, b \ge 0$,
$a^2 + b^2 = 1$.

Overall: Steps 1–3 look correct to me (Claude). Open questions are about presentation and Remark 2.

---

## 1. Step 2: can it be simplified with a convex combination + triangle inequality?

**Idea:** write $\Phi(x) = a^2\Phi(x) + b^2\Phi(x)$ and apply the triangle inequality.

**Why it fails.** $T_n$ is a weighted *sum* of independent variables, so its cdf is a
*convolution* of $F_1$ and $F_2$, not a mixture. $\Pr(aZ_1 + bZ_2 \le x)$ cannot be written as
$a^2(\cdot) + b^2(\cdot)$. The target bound $a^2\varepsilon_1 + b^2\varepsilon_2$ is false:

- $Z_2 \sim N(0,1)$ exactly, $Z_1 \sim N(0.2, 1)$, $a = b = 1/\sqrt2$.
- $\varepsilon_1 = 2\Phi(0.1) - 1 \approx 0.080$, $\varepsilon_2 = 0$, so $a^2\varepsilon_1 \approx 0.040$.
- But $T_n \sim N(0.141, 1)$, whose sup-distance from $\Phi$ is $\approx 0.056 > 0.040$.

Intuition: a shift in $Z_1$ passes into $T_n$ scaled by $a$, not $a^2$, and $a > a^2$.

**A genuine simplification: symmetric two-swap (telescoping) argument.**

$$
|\Pr(aZ_1+bZ_2\le x)-\Phi(x)|
\le \underbrace{|\Pr(aZ_1+bZ_2\le x)-\Pr(aW_1+bZ_2\le x)|}_{\text{condition on } Z_2:\ \le \varepsilon_1}
+ \underbrace{|\Pr(aW_1+bZ_2\le x)-\Pr(aW_1+bW_2\le x)|}_{\text{condition on } W_1:\ \le \varepsilon_2}
$$

The second term is $|\mathbb E[F_2((x-aW_1)/b) - \Phi((x-aW_1)/b)]| \le \varepsilon_2$, mirroring
the first. This removes the monotone-function / layer-cake parenthesis in the note. Handle
$a = 0$ or $b = 0$ separately (immediate).

### Does the convex-combination idea work with characteristic functions instead?

Partly, but not through a convex combination.

- Independence gives $\mathbb E[e^{itT_n}] = \varphi_1(at)\varphi_2(bt)$, and $a^2 + b^2 = 1$ gives
  $\varphi(t) = \varphi(at)\varphi(bt)$.
- Add and subtract a cross term, use $|\varphi_k| \le 1$:
  $|\varphi_1(at)\varphi_2(bt) - \varphi(at)\varphi(bt)| \le |\varphi_1(at)-\varphi(at)| + |\varphi_2(bt)-\varphi(bt)|$.
  This is telescoping with weights $1, 1$ (what the current Step 3 in `derivations.tex` does).
- The convex combination only appears in the Gaussian *exponent*:
  $\log\varphi(at) + \log\varphi(bt) = a^2\log\varphi(t) + b^2\log\varphi(t)$. That is just
  $a^2 + b^2 = 1$; it does not transfer to the errors, since $\log\varphi_1(at) \ne a^2\log\varphi_1(t)$
  unless $Z_1$ is exactly Gaussian.
- The bound $|\mathrm{cf}_T(t) - \varphi(t)| \le a^2\sup_{|u|\le|t|}|\varphi_1(u)-\varphi(u)| + b^2(\cdots)$
  is false. Same example as above: at $t = 1$, LHS $= 0.0857$, RHS $= 0.0606$; fails on about a
  third of $t \in (0, 5]$ (checked numerically).

**Why the cdf route is still better.** The cf route also needs (i) uniform convergence of cfs on
compacts (TODO citation still open in `derivations.tex`), (ii) Lévy continuity along each
sequence $\omega_n$, and (iii) a near-worst-case-sequence argument, because Lévy has no rate.
That is the old Steps 3, 4 and 7. The cdf swap gives the explicit bound $\varepsilon_1 + \varepsilon_2$
directly. Turning a cf bound into a cdf bound with a rate would need Esseen's smoothing inequality.

---

## 2. Theorems that could one-shot Step 2

**Multivariate Pólya (Ranga Rao, 1962).** If $P_n \Rightarrow P$ on $\mathbb R^d$ and $P$ is
absolutely continuous, then $\sup_C |P_n(C) - P(C)| \to 0$ over all measurable convex sets $C$.
Pólya is the $d = 1$ case. More general "uniformity class" version: Billingsley & Topsøe (1967).

Step 2 then reads:

1. $(Z_1, Z_2) \Rightarrow N(0, I_2)$ (independent components, joint cf factorises).
2. $\{T_n(\omega) \le x\} = \{(Z_1, Z_2) \in H\}$ with $H = \{u : au_1 + bu_2 \le x\}$ a half-plane (convex).
3. Ranga Rao gives $\sup_{\omega, x}|\Pr(T_n(\omega)\le x) - \Phi(x)| \to 0$. Uniformity in $\omega$ is free.

**Textbook-only alternative (compactness + subsequences).** Suppose the sup does not tend to 0;
take $\omega_n$ along a subsequence with sup $\ge \epsilon$; $(a_n, b_n)$ lie on the compact quarter
circle, so a further subsequence converges to $(a, b)$; $a_nZ_1 + b_nZ_2 = aZ_1 + bZ_2 + o_P(1)$
(tightness), so CMT + Slutsky give $\Rightarrow N(0,1)$; Pólya gives $\sup_x \to 0$, contradiction.

**Not useful:** Berry–Esseen needs $\hat L_k$ to be a normalised sum with finite third moments;
(A2) only assumes asymptotic normality.

**Trade-off.** Both shortcuts prove the error tends to 0 but lose the explicit bound
$\varepsilon_1 + \varepsilon_2$, which Remark 2 relies on.

**References to verify before citing (not in Zotero; theorem numbers and pages unchecked):**

- R. Ranga Rao (1962), "Relations between weak and uniform convergence of measures with
  applications", *Annals of Mathematical Statistics* 33, 659–680.
- P. Billingsley and F. Topsøe (1967), "Uniformity in weak convergence",
  *Z. Wahrscheinlichkeitstheorie verw. Gebiete* 7, 1–16.

---

## 3. Step 3: why not just apply Slutsky?

Slutsky at a fixed $\omega$ is not enough: the statement has $\inf_\omega$ *inside* the
$\liminf$, so it needs uniformity in $\omega$, and Slutsky is pointwise.

**Slutsky along a sequence does work:**

1. Pick deterministic $\omega_n$ with $\Pr(\hat T_n(\omega_n) \le z) \le \inf_\omega \Pr(\hat T_n(\omega) \le z) + 1/n$.
2. Step 2 gives $T_n(\omega_n) \Rightarrow N(0,1)$.
3. The weighted-average display gives $\hat v_n(\omega_n)/v_n(\omega_n) \overset{p}{\to} 1$.
4. Slutsky gives $\hat T_n(\omega_n) \Rightarrow N(0,1)$, hence $\liminf \inf_\omega \ge 1 - \alpha/2$.

This is essentially the old Steps 6–7.

**What the $\eta$-argument buys.** It is the proof of Slutsky written out with explicit error
terms, giving
$\Pr(\hat T_n(\omega) \le z) \ge \Phi(z\sqrt{1-\eta}) - \varepsilon_1 - \varepsilon_2 - \Pr(\delta_n > \eta)$
for every $n$ and $\omega$. Remark 2 needs that explicit bound. The Slutsky route only gives
pointwise-in-$P$ validity.

**Decision to make:** keep the $\eta$-argument if Remark 2 stays; use Slutsky-along-a-sequence if not.

---

## 4. Is Remark 2 (uniformity over distributions) correct?

My concern was: "we don't have uniformity over distributions; the proof relies on a normal
distribution." Claude's view: the remark is **correct but conditional**, and it should be reworded
rather than deleted.

**Where normality enters.** The data are never assumed normal; $F_k$ is arbitrary. Normality
appears only as (i) $\Phi$, the common limit in (A2), and (ii) the auxiliary $W_1, W_2$ used to
evaluate $\Pr(aW_1 + bW_2 \le x) = \Phi(x)$. Neither depends on the DGP $P$. A DGP enters the bound
only through $\varepsilon_1(P)$, $\varepsilon_2(P)$ and $\Pr_P(\delta_n > \eta)$.

**Why it holds.** The bound above is exact for every $n$, $\omega$ and $P$; Steps 1–2 are
finite-sample statements. So if those three quantities tend to 0 uniformly over a class
$\mathcal P$, then $\liminf_n \inf_{P\in\mathcal P}\inf_\omega \Pr_P(\ell \ge \hat A(\omega)) \ge 1 - \alpha/2$.

**Non-normal example where it holds.** $\hat L_k$ a sample mean of i.i.d. data, $\hat\sigma_k^2$
the sample variance, third moments $\le M$, variance $\ge c > 0$: Berry–Esseen gives
$\sup_{\mathcal P}\varepsilon_k \le CM/(c^{3/2}\sqrt{n_k}) \to 0$; with bounded fourth moments,
Chebyshev gives $\sup_{\mathcal P}\Pr(\delta_n > \eta) \to 0$. Skewed or heavy-tailed data within
those bounds are fine.

**Legitimate criticisms:**

1. (A2)–(A3) are pointwise in $P$, so the proposition as stated only proves pointwise validity.
   Uniform validity needs extra assumptions the note does not state (moment bounds, $\sigma_k^2$
   bounded away from 0).
2. Those assumptions do real work: Bernoulli($p$) data with $p \to 0$ allowed in the class breaks
   uniform $\varepsilon_k \to 0$, and uniform coverage can fail.
3. "The usual standard (Imbens & Manski 2004; Stoye 2009)" oversells it: those papers *assume*
   uniform asymptotic normality explicitly.

**Suggested rewording:** "Under uniform versions of (A2) and (A3) over a class $\mathcal P$ (e.g.
bounded third moments and variances bounded away from zero), the bound in Step 3 gives coverage
that is uniform over $\mathcal P$."

---

## To think about

- [ ] Do we want Remark 2 (uniformity over $P$) in the paper? This decides Step 2 (explicit swap
      vs Ranga Rao) and Step 3 ($\eta$-argument vs Slutsky).
- [ ] If keeping it: state uniform versions of (A2)–(A3) explicitly.
- [ ] If adopting the note: switch Step 2 to the symmetric two-swap version.
- [ ] Verify and add Ranga Rao (1962) / Billingsley & Topsøe (1967) to Zotero if citing.
