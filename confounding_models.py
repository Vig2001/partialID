import numpy as np
import scipy.integrate as integrate
import scipy.stats as stats
import matplotlib.pyplot as plt

# In this script I use OOP to make my simulations more systematic and easy to follow

class BaseConfoundingModel:
    def __init__(self, alpha=1.0, gamma_uy=2.0, gamma_ut=4.0):
        self.alpha = alpha
        self.gamma_uy = gamma_uy
        self.gamma_ut = gamma_ut

    def u_mean(self, X): return 0.0
    def u_std(self, X): return 1.0

    def pdf_u_given_x(self, u, X):
        """Default to standard Gaussian. Override this in subclasses for different distributions."""
        return stats.norm.pdf(u, loc=self.u_mean(X), scale=self.u_std(X))

    def integration_bounds(self, X):
        """
        Default to +/- 10 std deviations for Gaussian. 
        Override this for bounded distributions (like Uniform or Beta).
        """
        mu = self.u_mean(X)
        sig = self.u_std(X)
        return mu - 10 * sig, mu + 10 * sig

    def potential_outcomes(self, X, U):
        raise NotImplementedError("Subclasses must define this")
        
    def propensity_score(self, X, U):
        raise NotImplementedError("Subclasses must define this")

    def marginal_propensity(self, X, T):
        lower, upper = self.integration_bounds(X)
        if T == 1:
            return integrate.quad(lambda u: self.propensity_score(X, u) * self.pdf_u_given_x(u, X), lower, upper)[0]
        else:
            return integrate.quad(lambda u: (1 - self.propensity_score(X, u)) * self.pdf_u_given_x(u, X), lower, upper)[0]

    def catt(self, X):
        lower, upper = self.integration_bounds(X)
        marg_p1 = self.marginal_propensity(X, T=1)
        return integrate.quad(
            lambda u: self.potential_outcomes(X, u)[2] * (self.propensity_score(X, u) / marg_p1) * self.pdf_u_given_x(u, X), 
            lower, upper
        )[0]

    def catc(self, X):
        lower, upper = self.integration_bounds(X)
        marg_p0 = self.marginal_propensity(X, T=0)
        return integrate.quad(
            lambda u: self.potential_outcomes(X, u)[2] * ((1 - self.propensity_score(X, u)) / marg_p0) * self.pdf_u_given_x(u, X), 
            lower, upper
        )[0]

    def true_cate(self, X):
        cate_1 = self.catt(X) * self.marginal_propensity(X, T=1)
        cate_0 = self.catc(X) * self.marginal_propensity(X, T=0)
        return cate_1 + cate_0

    def os_cate(self, X):
        lower, upper = self.integration_bounds(X)
        marg_p1 = self.marginal_propensity(X, T=1)
        marg_p0 = self.marginal_propensity(X, T=0)
        
        e_y1_given_t1 = integrate.quad(
            lambda u: self.potential_outcomes(X, u)[1] * (self.propensity_score(X, u) / marg_p1) * self.pdf_u_given_x(u, X), 
            lower, upper
        )[0]
        
        e_y0_given_t0 = integrate.quad(
            lambda u: self.potential_outcomes(X, u)[0] * ((1 - self.propensity_score(X, u)) / marg_p0) * self.pdf_u_given_x(u, X), 
            lower, upper
        )[0]
        
        return e_y1_given_t1 - e_y0_given_t0
    
# Strictly Linear
class LinearModel(BaseConfoundingModel):
    def potential_outcomes(self, X, U):
        Y0 = 1.0 + 0.5 * X + U
        Y1 = 1.0 + 2.0 * X + self.gamma_uy * U
        return Y0, Y1, Y1 - Y0
        
    def propensity_score(self, X, U):
        logit = self.alpha * X + self.gamma_ut * U
        return 1 / (1 + np.exp(-logit))

# Linear combinations of Nonlinear terms
class NonlinearFeaturesModel(BaseConfoundingModel):
    def potential_outcomes(self, X, U):
        Y0 = 1.0 + 0.5 * (X**2) + U
        Y1 = 1.0 + 2.0 * (X**2) + self.gamma_uy * (U**2)
        return Y0, Y1, Y1 - Y0
        
    def propensity_score(self, X, U):
        logit = self.alpha * (X**2) + self.gamma_ut * np.sin(U)
        return 1 / (1 + np.exp(-logit))

# Completely Nonlinear / Interacting
class ComplexNonlinearModel(BaseConfoundingModel):
    def u_mean(self, X): 
        return 0.0 if X < 1 else 4.0

    def potential_outcomes(self, X, U):
        Y0 = np.exp(0.5 * X) + U
        Y1 = np.exp(X) + self.gamma_uy * X * U  # Interaction between X and U
        return Y0, Y1, Y1 - Y0
        
    def propensity_score(self, X, U):
        logit = self.alpha * X + self.gamma_ut * (X * U) # Interaction in treatment assignment
        return 1 / (1 + np.exp(-logit))