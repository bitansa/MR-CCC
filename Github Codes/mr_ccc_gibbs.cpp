// =============================================================================
// mr_ccc_gibbs.cpp
//
// Gibbs sampler for MR-CCC: Bayesian Mendelian Randomization for Causal
// Cell-Cell Communication.
//
// -----------------------------------------------------------------------------
// MODEL OVERVIEW (single ordered pair, scalar X / Z / Y per donor)
// -----------------------------------------------------------------------------
//
//  First stage - sender:
//    X = G * Pi_X  +  V * Alpha_X  +  eps_X,    eps_X ~ N(0, sigma_X^2 * I)
//
//  First stage - receiver:
//    Z = H * Pi_Z  +  V * Alpha_Z  +  eps_Z,    eps_Z ~ N(0, sigma_Z^2 * I)
//
//  Second stage - outcome:
//    Y = mu  +  X* * Beta_X  +  X*Z* * Beta_XZ
//           +  Z* * Beta_Z  +  V * Alpha_Y  +  eps_Y,
//                                               eps_Y ~ N(0, sigma_Y^2 * I)
//
//  where X* = G*Pi_X + V*Alpha_X  and  Z* = H*Pi_Z + V*Alpha_Z
//  are the IV-projected (genetically predicted) expression levels.
//
//  Spike-and-slab prior on the causal block (Beta_X, Beta_XZ), indexed by
//  the inclusion indicator gamma:
//    gamma = 1 (slab):  (Beta_X, Beta_XZ) ~ N(0, gBeta * sigma_Y^2 * D0^{-1})
//    gamma = 0 (spike): same form but variance scaled down by nu1  (nu1 << 1)
//  and on the receptor main effect:
//    Beta_Z ~ N(0, gZ * sigma_Y^2 / d_Z)
//
//  D0 = diag(||X^*||^2, ||X^* o Z^*||^2) and d_Z = ||Z^*||^2 are computed ONCE
//  from the plug-in first stage (least squares of X on [G V] and of Z on
//  [H V]) and held fixed. They scale the prior by the size of each regressor,
//  so the prior is placed on the size of the effect on Y rather than on raw
//  coefficients whose scale is set by units and by instrument strength.
//  Equivalently: independent N(0, gBeta * s_gamma * sigma_Y^2) priors on the
//  coefficients of the standardised regressors X^*/||X^*|| and
//  X^*Z^*/||X^*Z^*||, with the standardisation fixed at the plug-in fit.
//
//  The scale is fixed rather than recomputed from the current draws because
//  X^* and Z^* are functions of the first-stage parameters. A prior whose
//  covariance depends on those parameters would add a determinant and a
//  quadratic-form term to their full conditionals; the first-stage updates
//  below do not carry such terms, and with a fixed scale none is needed.
//  The plug-in values are a scale, not an estimate: X^* and Z^* in the
//  likelihood are still rebuilt from the current draws at every sweep.
//
//  gamma therefore acts through the PRIOR VARIANCE of the causal block and
//  does not multiply the mean function: under the spike the two coefficients
//  are drawn from a distribution concentrated tightly at zero rather than set
//  to zero exactly. This is the continuous spike-and-slab formulation; with
//  nu1 = 1e-4 it is numerically indistinguishable from a point-mass spike,
//  and P(gamma = 1 | data) is the posterior inclusion probability.
//
// -----------------------------------------------------------------------------
// INPUTS
// -----------------------------------------------------------------------------
//  X       (n x 1)   Ligand expression in the sender cell type
//  Z       (n x 1)   Receptor expression in the receiver cell type
//  Y       (n x 1)   Downstream pathway activity in the receiver cell type
//  G       (n x pG)  Sender cis-eQTL genotypes (instruments for X)
//  H       (n x pH)  Receiver cis-eQTL genotypes (instruments for Z)
//  V       (n x pV)  Shared covariates (e.g., ancestry PCs, cell proportions)
//
// MCMC SETTINGS
//  n_iter  Total number of Gibbs iterations
//  burn_in Burn-in length (iterations discarded before accumulation)
//  thin    Thinning interval (retain every `thin`-th post-burn-in sample)
//
// HYPERPARAMETERS
//  a_sigma, b_sigma  InvGamma shape / scale for sigma_X^2, sigma_Z^2, sigma_Y^2
//  a_rho,   b_rho    Beta shape parameters for rho (prior inclusion probability)
//  nu1               Spike variance multiplier; nu1 << 1 enforces near-zero spike
//  gG, gH, gV, gZ, gBeta   g-prior scale factors for each coefficient block
//  ridge             Small diagonal perturbation for numerical stability
//
// OUTPUT
//  Named R list containing:
//
//   (a) Posterior MEANS for all model parameters, including
//       gamma_mean = P(gamma = 1 | data), the posterior inclusion
//       probability (PIP) used as the communication score.
//
//   (b) Retained DRAWS for the scalar parameters of inferential interest:
//         n_keep          number of retained draws
//         Beta_X_draws    ligand main effect
//         Beta_XZ_draws   receptor-modulated interaction effect
//         Beta_Z_draws    receptor main effect
//         gamma_draws     inclusion indicator (0/1) per retained iteration
//         mu_draws        intercept
//         loglik_draws    joint log-likelihood of the three equations at the
//                         retained draw (likelihood only; no prior terms)
//
//   (c) prior_scale = c(d_X, d_XZ, d_Z), the fixed plug-in scale used by the
//       priors on (Beta_X, Beta_XZ) and Beta_Z, so that every fit records the
//       scale it was run under.
//
//       Draws are required for quantities that cannot be recovered from
//       posterior means: credible intervals for functionals such as the
//       sign-reversal threshold tau = -Beta_X / Beta_XZ (a ratio, so it
//       needs paired draws), Monte Carlo standard errors and effective
//       sample size for the PIP, and trace plots for convergence checks.
// =============================================================================

// -----------------------------------------------------------------------------
// PROGRESS REPORTING: verbose
//  When verbose is true, a progress line is written to the R console every
//  1000 iterations. The default is silent, because a full analysis runs
//  several chains on each of many triplets.
//
// OPTIONAL INITIALISATION CONTROLS: init_gamma, init_scale
// -----------------------------------------------------------------------------
//  The Gelman-Rubin diagnostic compares chains started from OVERDISPERSED
//  values. Chains that differ only in their random number stream share a
//  starting point and therefore provide a weaker check on convergence. These
//  two arguments allow dispersed starts:
//
//    init_gamma   starting value of the inclusion indicator (0 = spike,
//                 1 = slab). Running some chains from each is the natural
//                 dispersion for a spike-and-slab model.
//    init_scale   if greater than zero, Beta_X and Beta_XZ are initialised
//                 as draws from a normal with mean 0 and STANDARD DEVIATION
//                 init_scale, instead of at zero. (R::rnorm takes the
//                 standard deviation, not the variance.)
//
//  The defaults (init_gamma = 1, init_scale = 0) give a deterministic start
//  with the slab active and both causal coefficients at zero.
//
//  NOTE ON THE ARGUMENT LIST: comments inside the argument list are kept to
//  a single trailing comment per parameter. Rcpp flattens the signature onto
//  one line when generating its wrapper, and a multi-line comment block placed
//  among the parameters would be merged into one line and comment out the
//  remainder, producing a compile error.
// -----------------------------------------------------------------------------

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// =============================================================================
// Helper functions
// =============================================================================

// Draw one sample from InvGamma(a, b).
// Uses the identity: if X ~ Gamma(a, 1/b) then 1/X ~ InvGamma(a, b).
inline double rinvgamma_1(double a, double b) {
  return 1.0 / R::rgamma(a, 1.0 / b);
}

// Numerically stable log-sum-exp for exactly two terms:
//   log( exp(a) + exp(b) ) = m + log( exp(a-m) + exp(b-m) ),  m = max(a, b)
inline double logsumexp2(double a, double b) {
  double m = (a > b) ? a : b;
  return m + std::log(std::exp(a - m) + std::exp(b - m));
}

// Inverse of a symmetric positive-definite matrix, stopping with an
// informative message instead of a raw Armadillo error when the matrix is
// numerically singular (for example, collinear instruments or covariates).
inline arma::mat safe_inv_sympd(const arma::mat& A, const char* what) {
  arma::mat out;
  if (!arma::inv_sympd(out, A)) {
    Rcpp::stop("Numerical failure inverting the %s matrix; check for "
               "collinear instruments or covariates, or increase 'ridge'.",
               what);
  }
  return out;
}


// =============================================================================
// Main exported function
// =============================================================================
//
// The roxygen block below (lines beginning //') is copied verbatim into
// R/RcppExports.R by Rcpp::compileAttributes(), which is how the package
// documentation for this function is generated. It lives here rather than in
// RcppExports.R so that regenerating the wrapper cannot delete it.

//' Gibbs sampler for the MR-CCC model
//'
//' Runs the blocked Gibbs sampler for a single ligand-receptor-pathway
//' triplet and returns posterior means together with the retained draws of
//' the scalar parameters of inferential interest. This is the low-level
//' engine called by [mr_ccc()]; most users should call [mr_ccc()], which
//' validates and centres the inputs, runs several chains, and computes
//' convergence diagnostics.
//'
//' @param X Numeric matrix of dimension n x 1: centred ligand expression in
//'   the sender cell type.
//' @param Z Numeric matrix of dimension n x 1: centred receptor expression in
//'   the receiver cell type.
//' @param Y Numeric matrix of dimension n x 1: centred pathway activity in
//'   the receiver cell type.
//' @param G Numeric matrix of dimension n x pG: column-centred sender
//'   cis-eQTL genotypes (instruments for X).
//' @param H Numeric matrix of dimension n x pH: column-centred receiver
//'   cis-eQTL genotypes (instruments for Z).
//' @param V Numeric matrix of dimension n x pV, pV >= 1: column-centred
//'   shared covariates. The first stages have no intercept, so all three
//'   design matrices must be centred, as [mr_ccc()] does.
//' @param n_iter Total number of Gibbs iterations.
//' @param burn_in Number of initial iterations discarded.
//' @param thin Thinning interval; every `thin`-th post-burn-in iteration is
//'   retained.
//' @param a_sigma,b_sigma Inverse-gamma shape and scale for the three
//'   residual variances.
//' @param a_rho,b_rho Beta shape parameters for the prior inclusion
//'   probability rho.
//' @param nu1 Spike variance multiplier; small values enforce a near-zero
//'   spike.
//' @param gG,gV,gH,gZ,gBeta Prior scale factors for the first-stage sender
//'   effects, the covariate effects, the first-stage receiver effects, the
//'   receptor main effect, and the causal block (Beta_X, Beta_XZ). The first
//'   three are Zellner g-priors on the observed designs G, V and H; the last
//'   two scale a fixed diagonal computed once from the least-squares first
//'   stage (see Details).
//' @param ridge Diagonal ridge added before every matrix inversion.
//' @param init_gamma Starting value (0 or 1) of the inclusion indicator.
//' @param init_scale If greater than zero, Beta_X and Beta_XZ start from
//'   independent normal draws with standard deviation `init_scale`.
//' @param verbose Logical; if `TRUE`, a progress line is printed every 1000
//'   iterations.
//'
//' @details
//' The priors on the second-stage coefficients are
//' \deqn{(\beta_X, \beta_{XZ}) \mid \gamma, \sigma_Y^2 \sim
//'   N_2(0,\; g_\beta s_\gamma \sigma_Y^2 D_0^{-1}), \qquad
//'   \beta_Z \mid \sigma_Y^2 \sim N(0,\; g_Z \sigma_Y^2 / d_Z),}
//' with \eqn{D_0 = \mathrm{diag}(\|\hat X^*\|^2, \|\hat X^* \circ \hat Z^*\|^2)}
//' and \eqn{d_Z = \|\hat Z^*\|^2} computed once from least squares of X on
//' \[G V\] and of Z on \[H V\] and held fixed, and \eqn{s_\gamma = 1} in the
//' slab and `nu1` in the spike. Because the projected regressors \eqn{X^*}
//' and \eqn{Z^*} are functions of the first-stage parameters, a prior whose
//' covariance was recomputed from the current draws would enter the full
//' conditionals of those parameters; fixing the scale keeps every update an
//' exact conjugate step. The values used are returned in `prior_scale`.
//'
//' @return A named list with posterior means (`Pi_X_mean`, `Alpha_X_mean`,
//'   `Pi_Z_mean`, `Alpha_Z_mean`, `Beta_X_mean`, `Beta_XZ_mean`,
//'   `Beta_Z_mean`, `Alpha_Y_mean`, `mu_mean`, `sigma_X_sq_mean`,
//'   `sigma_Z_sq_mean`, `sigma_Y_sq_mean`, `gamma_mean`, `rho_mean`), the
//'   number of retained draws `n_keep`, the retained draws `Beta_X_draws`,
//'   `Beta_XZ_draws`, `Beta_Z_draws`, `gamma_draws`, `mu_draws` and
//'   `loglik_draws` (each a numeric vector of length `n_keep`), and
//'   `prior_scale`, a named numeric vector giving `d_X`, `d_XZ` and `d_Z`.
//'   `loglik_draws` is the joint
//'   log-likelihood of the three equations evaluated at each retained draw,
//'   excluding prior terms.
//'
//' @seealso [mr_ccc()] for the user-facing interface.
//'
//' @examples
//' sim <- simulate_mrccc(n = 120, seed = 1)
//' ctr <- function(v) matrix(v - mean(v), ncol = 1)
//' ctr_cols <- function(M) scale(M, center = TRUE, scale = FALSE)
//' out <- mr_ccc_gibbs(ctr(sim$X), ctr(sim$Z), ctr(sim$Y),
//'                     ctr_cols(sim$G), ctr_cols(sim$H), ctr_cols(sim$V),
//'                     n_iter = 400, burn_in = 100, thin = 1)
//' out$gamma_mean
//' out$prior_scale
//' length(out$Beta_XZ_draws) == out$n_keep
//'
//' @export
// [[Rcpp::export]]
List mr_ccc_gibbs(
    const arma::mat& X,       // (n x 1) ligand expression - sender
    const arma::mat& Z,       // (n x 1) receptor expression - receiver
    const arma::mat& Y,       // (n x 1) pathway activity - receiver
    const arma::mat& G,       // (n x pG) sender cis-eQTL instruments
    const arma::mat& H,       // (n x pH) receiver cis-eQTL instruments
    const arma::mat& V,       // (n x pV) shared covariates
    int    n_iter   = 20000,  // total MCMC iterations
    int    burn_in  = 2000,   // burn-in (discarded)
    int    thin     = 1,      // thinning interval (1 = keep every draw)
    double a_sigma  = 3.0,    // InvGamma shape for variance parameters
    double b_sigma  = 2.0,    // InvGamma scale for variance parameters
    double a_rho    = 3.0,    // Beta shape (a) for inclusion probability rho
    double b_rho    = 1.0,    // Beta shape (b) for inclusion probability rho
    double nu1      = 1e-4,   // spike variance multiplier
    double gG       = 100.0,  // g-prior scale for Pi_X  (sender eQTL effects)
    double gV       = 100.0,  // g-prior scale for Alpha_X, Alpha_Z, Alpha_Y
    double gH       = 100.0,  // g-prior scale for Pi_Z  (receiver eQTL effects)
    double gZ       = 100.0,  // g-prior scale for Beta_Z (receptor main effect)
    double gBeta    = 100.0,  // g-prior scale for (Beta_X, Beta_XZ) (causal block)
    double ridge    = 1e-8,   // diagonal ridge for stable matrix inversion
    int    init_gamma = 1,    // starting value of the inclusion indicator
    double init_scale = 0.0,  // if > 0, Beta_X, Beta_XZ start ~ N(0, sd=init_scale)
    bool   verbose = false    // if true, report progress every 1000 iterations
) {

  // --------------------------------------------------------------------------
  // Dimensions
  // --------------------------------------------------------------------------
  const int n  = X.n_rows;
  const int pG = G.n_cols;
  const int pH = H.n_cols;
  const int pV = V.n_cols;

  // --------------------------------------------------------------------------
  // Input checks. mr_ccc() validates its inputs before calling this function;
  // these checks protect direct calls, where a bad value would otherwise end
  // in an obscure linear-algebra error or, for thin = 0, a division by zero.
  // --------------------------------------------------------------------------
  if (n < 2 || X.n_cols != 1 || Z.n_cols != 1 || Y.n_cols != 1)
    Rcpp::stop("'X', 'Z' and 'Y' must be n x 1 matrices with n >= 2.");
  if ((int) Z.n_rows != n || (int) Y.n_rows != n || (int) G.n_rows != n ||
      (int) H.n_rows != n || (int) V.n_rows != n)
    Rcpp::stop("All inputs must have the same number of rows.");
  if (pG < 1 || pH < 1 || pV < 1)
    Rcpp::stop("'G', 'H' and 'V' must each have at least one column.");
  if (n_iter < 1 || burn_in < 0 || burn_in >= n_iter)
    Rcpp::stop("Require n_iter >= 1 and 0 <= burn_in < n_iter.");
  if (thin < 1)
    Rcpp::stop("'thin' must be an integer >= 1.");
  if (!(nu1 > 0.0 && nu1 < 1.0))
    Rcpp::stop("'nu1' must lie strictly between 0 and 1.");
  if (!(a_sigma > 0.0 && b_sigma > 0.0 && a_rho > 0.0 && b_rho > 0.0))
    Rcpp::stop("'a_sigma', 'b_sigma', 'a_rho' and 'b_rho' must be positive.");
  if (!(gG > 0.0 && gV > 0.0 && gH > 0.0 && gZ > 0.0 && gBeta > 0.0) ||
      ridge < 0.0)
    Rcpp::stop("The g scales must be positive and 'ridge' non-negative.");

  // --------------------------------------------------------------------------
  // Precompute Gram matrices and regularized inverses.
  // These are fixed throughout the MCMC and reused at every iteration.
  // --------------------------------------------------------------------------
  const arma::mat G_Star = G.t() * G;  // (pG x pG)
  const arma::mat V_Star = V.t() * V;  // (pV x pV)
  const arma::mat H_Star = H.t() * H;  // (pH x pH)

  const arma::mat G_Star_Inv = safe_inv_sympd(G_Star + ridge * arma::eye(pG, pG), "G'G");
  const arma::mat V_Star_Inv = safe_inv_sympd(V_Star + ridge * arma::eye(pV, pV), "V'V");
  const arma::mat H_Star_Inv = safe_inv_sympd(H_Star + ridge * arma::eye(pH, pH), "H'H");

  // Projection matrices (X'X)^{-1} X' -- appear in g-prior posterior means
  const arma::mat G_Tilde = G_Star_Inv * G.t();  // (pG x n)
  const arma::mat V_Tilde = V_Star_Inv * V.t();  // (pV x n)
  const arma::mat H_Tilde = H_Star_Inv * H.t();  // (pH x n)

  // --------------------------------------------------------------------------
  // Initialize model parameters
  // --------------------------------------------------------------------------

  // First-stage sender:   X = G * Pi_X + V * Alpha_X + eps_X
  arma::rowvec Pi_X    = arma::zeros<arma::rowvec>(pG);  // eQTL effects on ligand
  arma::rowvec Alpha_X = arma::zeros<arma::rowvec>(pV);  // covariate effects on ligand

  // First-stage receiver: Z = H * Pi_Z + V * Alpha_Z + eps_Z
  arma::rowvec Pi_Z    = arma::zeros<arma::rowvec>(pH);  // eQTL effects on receptor
  arma::rowvec Alpha_Z = arma::zeros<arma::rowvec>(pV);  // covariate effects on receptor

  // Second-stage causal effects
  double Beta_X  = 0.0;  // main causal effect:  ligand X* -> pathway Y
  double Beta_XZ = 0.0;  // interaction effect:  X* modulated by Z* -> Y
  double Beta_Z  = 0.0;  // direct receptor effect: Z* -> Y (non-causal pathway)

  // Joint row-vector for (Beta_X, Beta_XZ) -- used in spike-and-slab and sigma_Y updates
  arma::rowvec Beta_row(2);
  Beta_row(0) = Beta_X;
  Beta_row(1) = Beta_XZ;

  // Covariate effects on the outcome
  arma::rowvec Alpha_Y = arma::zeros<arma::rowvec>(pV);

  // Intercept for Y.  Prior: mu ~ N(0, c_mu * sigma_Y^2),  c_mu = n (weakly informative).
  double mu   = 0.0;
  double c_mu = static_cast<double>(n);

  // Variance parameters -- initialized to (biased) sample variances
  double sigma_X_sq = arma::as_scalar(arma::var(X.col(0))) * (n - 1.0) / n;
  double sigma_Z_sq = arma::as_scalar(arma::var(Z.col(0))) * (n - 1.0) / n;
  double sigma_Y_sq = arma::as_scalar(arma::var(Y.col(0))) * (n - 1.0) / n;

  // Spike-and-slab inclusion indicator and prior inclusion probability
  int    gamma = 1;    // gamma = 1: slab (active signal); gamma = 0: spike (null)
  double rho   = 0.5;  // prior inclusion probability, updated each iteration

  // ---- Apply optional dispersed initialisation ------------------------------
  // With the defaults (init_gamma = 1, init_scale = 0) the chain starts with the
  // slab active and both causal coefficients at zero. Other values give the
  // overdispersed starts required for a meaningful Gelman-Rubin statistic.
  gamma = (init_gamma == 0) ? 0 : 1;
  if (init_scale > 0.0) {
    Beta_X  = R::rnorm(0.0, init_scale);
    Beta_XZ = R::rnorm(0.0, init_scale);
  }

  // --------------------------------------------------------------------------
  // Accumulators for posterior means (post-burn-in, thinned samples only)
  // --------------------------------------------------------------------------
  arma::rowvec Pi_X_sum    = arma::zeros<arma::rowvec>(pG);
  arma::rowvec Alpha_X_sum = arma::zeros<arma::rowvec>(pV);
  arma::rowvec Pi_Z_sum    = arma::zeros<arma::rowvec>(pH);
  arma::rowvec Alpha_Z_sum = arma::zeros<arma::rowvec>(pV);
  arma::rowvec Alpha_Y_sum = arma::zeros<arma::rowvec>(pV);

  double Beta_X_sum  = 0.0, Beta_XZ_sum = 0.0, Beta_Z_sum = 0.0;
  double sX_sum      = 0.0, sZ_sum      = 0.0, sY_sum     = 0.0;
  double gamma_sum   = 0.0, rho_sum     = 0.0, mu_sum     = 0.0;

  // --------------------------------------------------------------------------
  // Retained draws for the scalar parameters of inferential interest.
  //
  // Posterior MEANS alone are not sufficient for several quantities:
  //
  //   * the sign-reversal threshold  tau = -Beta_X / Beta_XZ  is a ratio of two
  //     parameters, so its posterior requires PAIRED draws, not two means;
  //   * Monte Carlo standard error and effective sample size for the PIP
  //     require the gamma chain itself;
  //   * trace plots and convergence diagnostics require the chains.
  //
  // Only the six scalars needed for these purposes are stored, so the memory
  // cost is negligible (n_keep x 6 doubles; < 10 MB even for 200,000 kept
  // draws). The high-dimensional blocks (Pi_X, Pi_Z, Alpha_*) are still
  // summarised by their means only.
  //
  // The sixth scalar is the joint log-likelihood of the three equations,
  // evaluated at each retained draw. It is a single scalar summary of overall
  // model fit and is the conventional quantity to display in a trace plot when
  // the parameter space is too large to plot component by component.
  // --------------------------------------------------------------------------
  const int n_keep = (n_iter > burn_in && thin > 0)
                       ? ((n_iter - burn_in) / thin) : 0;

  arma::vec Beta_X_draws  = arma::zeros<arma::vec>(n_keep);
  arma::vec Beta_XZ_draws = arma::zeros<arma::vec>(n_keep);
  arma::vec Beta_Z_draws  = arma::zeros<arma::vec>(n_keep);
  arma::vec gamma_draws   = arma::zeros<arma::vec>(n_keep);
  arma::vec mu_draws      = arma::zeros<arma::vec>(n_keep);
  arma::vec loglik_draws  = arma::zeros<arma::vec>(n_keep);

  // Constant term of a univariate Gaussian log-density, used below.
  const double LOG_2PI = std::log(2.0 * M_PI);

  // --------------------------------------------------------------------------
  // Column vectors extracted once (avoids repeated .col(0) inside the loop)
  // --------------------------------------------------------------------------
  const arma::vec x = X.col(0);
  const arma::vec z = Z.col(0);
  const arma::vec y = Y.col(0);

  // --------------------------------------------------------------------------
  // Fixed prior scale for the second-stage coefficient blocks.
  //
  // Plug-in first stage: least squares of X on [G V] and of Z on [H V],
  // regularised with the same ridge used for the Gram matrices above. The
  // squared norms of the resulting projected regressors set the prior
  // variance of (Beta_X, Beta_XZ) and Beta_Z for the whole run.
  //
  // These are computed from the observed data alone and never revisited, so
  // the priors they scale do not depend on any parameter the sampler updates.
  // That is what makes the first-stage updates in Steps 1, 2, 4 and 5 exact:
  // a prior whose covariance moved with the current X*, Z* would add a
  // determinant and a quadratic-form term to those conditionals, and the
  // updates below carry neither. The model overview at the top of this file
  // gives the equivalent statement in terms of standardised regressors.
  //
  // The ridge added to each squared norm keeps the prior proper if a plug-in
  // regressor is numerically zero, which can happen for an interaction column
  // when one first stage has no signal at all.
  // --------------------------------------------------------------------------
  const arma::mat WG = arma::join_horiz(G, V);   // (n x (pG + pV))
  const arma::mat WH = arma::join_horiz(H, V);   // (n x (pH + pV))
  const arma::vec theta_X = arma::solve(
    WG.t() * WG + ridge * arma::eye(pG + pV, pG + pV), WG.t() * x);
  const arma::vec theta_Z = arma::solve(
    WH.t() * WH + ridge * arma::eye(pH + pV, pH + pV), WH.t() * z);
  const arma::vec X_hat  = WG * theta_X;   // plug-in X*
  const arma::vec Z_hat  = WH * theta_Z;   // plug-in Z*
  const arma::vec XZ_hat = X_hat % Z_hat;  // plug-in X* o Z*

  const double d_X  = arma::dot(X_hat,  X_hat)  + ridge;
  const double d_XZ = arma::dot(XZ_hat, XZ_hat) + ridge;
  const double d_Z  = arma::dot(Z_hat,  Z_hat)  + ridge;

  arma::mat D0 = arma::zeros<arma::mat>(2, 2);   // prior precision shape for (Beta_X, Beta_XZ)
  D0(0, 0) = d_X;
  D0(1, 1) = d_XZ;

  // ==========================================================================
  // MCMC (Gibbs) loop
  // ==========================================================================
  int save_idx = 0;

  for (int it = 1; it <= n_iter; ++it) {

    // Allow long runs to be interrupted from R.
    if (it % 1000 == 0) Rcpp::checkUserInterrupt();

    // IV-projected receiver expression and the resulting coefficient of X in Y
    // (computed from parameters at the start of the iteration)
    arma::vec Z_Dash  = H * Pi_Z.t() + V * Alpha_Z.t();  // Z* = IV-projected receptor  (n x 1)
    arma::vec WX_diag = Beta_X + Beta_XZ * Z_Dash;        // dY/dX* at each obs          (n x 1)

    // =========================================================================
    // BLOCK 1 -- Sender first stage: Pi_X, Alpha_X, sigma_X^2
    // =========================================================================

    // ---- Step 1: Pi_X | rest ------------------------------------------------
    // Posterior is MVN. The precision combines the X-equation g-prior and the
    // Y-equation likelihood (X* = G*Pi_X + V*Alpha_X enters Y through WX_diag).
    arma::mat GW = G.each_col() % WX_diag;  // G weighted by dY/dX*  (n x pG)

    arma::mat Sigma_Pi_X = safe_inv_sympd(
      (1.0 / sigma_X_sq) * (1.0 + 1.0 / gG) * G_Star +
      (1.0 / sigma_Y_sq) * (GW.t() * GW) +
      ridge * arma::eye(pG, pG), "Pi_X precision");

    arma::vec rhs_PiX =
      (1.0 / sigma_X_sq) * G.t() * (x - V * Alpha_X.t()) +
      (1.0 / sigma_Y_sq) * GW.t() * (y - mu
        - Z_Dash * Beta_Z
        - V * Alpha_Y.t()
        - (V * Alpha_X.t()) % WX_diag);

    arma::rowvec m_Pi_X = (Sigma_Pi_X * rhs_PiX).t();
    Pi_X = arma::mvnrnd(m_Pi_X.t(), Sigma_Pi_X, 1).t();

    // ---- Step 2: Alpha_X | rest ---------------------------------------------
    // Same structure as Pi_X but for the covariate block V.
    arma::mat VW = V.each_col() % WX_diag;  // V weighted by dY/dX*  (n x pV)

    arma::mat Sigma_Alpha_X = safe_inv_sympd(
      (1.0 / sigma_X_sq) * (1.0 + 1.0 / gV) * V_Star +
      (1.0 / sigma_Y_sq) * (VW.t() * VW) +
      ridge * arma::eye(pV, pV), "Alpha_X precision");

    arma::vec rhs_AlphaX =
      (1.0 / sigma_X_sq) * V.t() * (x - G * Pi_X.t()) +
      (1.0 / sigma_Y_sq) * VW.t() * (y - mu
        - Z_Dash * Beta_Z
        - V * Alpha_Y.t()
        - (G * Pi_X.t()) % WX_diag);

    arma::rowvec m_Alpha_X = (Sigma_Alpha_X * rhs_AlphaX).t();
    Alpha_X = arma::mvnrnd(m_Alpha_X.t(), Sigma_Alpha_X, 1).t();

    // ---- Step 3: sigma_X^2 | rest -------------------------------------------
    // Conjugate InvGamma. Scale includes the g-prior quadratic penalty on Pi_X
    // and Alpha_X so that the marginal prior on X integrates correctly.
    arma::vec x_res = x - (G * Pi_X.t() + V * Alpha_X.t());
    double penX = arma::as_scalar(Pi_X * G_Star * Pi_X.t()) / gG
      + arma::as_scalar(Alpha_X * V_Star * Alpha_X.t()) / gV;
    double sX_shape = a_sigma + 0.5 * (n + pG + pV);
    double sX_scale = b_sigma + 0.5 * (arma::dot(x_res, x_res) + penX);
    sigma_X_sq = rinvgamma_1(sX_shape, sX_scale);

    // Updated IV-projected sender expression (used in Z block and Y block)
    arma::vec X_Dash  = G * Pi_X.t() + V * Alpha_X.t();  // X* = IV-projected ligand  (n x 1)
    arma::vec WZ_diag = Beta_Z + Beta_XZ * X_Dash;         // dY/dZ* at each obs        (n x 1)

    // =========================================================================
    // BLOCK 2 -- Receiver first stage: Pi_Z, Alpha_Z, sigma_Z^2
    // =========================================================================

    // ---- Step 4: Pi_Z | rest ------------------------------------------------
    arma::mat HW = H.each_col() % WZ_diag;  // H weighted by dY/dZ*  (n x pH)

    arma::mat Sigma_Pi_Z = safe_inv_sympd(
      (1.0 / sigma_Z_sq) * (1.0 + 1.0 / gH) * H_Star +
      (1.0 / sigma_Y_sq) * (HW.t() * HW) +
      ridge * arma::eye(pH, pH), "Pi_Z precision");

    arma::vec rhs_PiZ =
      (1.0 / sigma_Z_sq) * H.t() * (z - V * Alpha_Z.t()) +
      (1.0 / sigma_Y_sq) * HW.t() * (y - mu
        - X_Dash * Beta_X
        - V * Alpha_Y.t()
        - (V * Alpha_Z.t()) % WZ_diag);

    arma::rowvec m_Pi_Z = (Sigma_Pi_Z * rhs_PiZ).t();
    Pi_Z = arma::mvnrnd(m_Pi_Z.t(), Sigma_Pi_Z, 1).t();

    // ---- Step 5: Alpha_Z | rest ---------------------------------------------
    arma::mat VWz = V.each_col() % WZ_diag;  // V weighted by dY/dZ*  (n x pV)

    arma::mat Sigma_Alpha_Z = safe_inv_sympd(
      (1.0 / sigma_Z_sq) * (1.0 + 1.0 / gV) * V_Star +
      (1.0 / sigma_Y_sq) * (VWz.t() * VWz) +
      ridge * arma::eye(pV, pV), "Alpha_Z precision");

    arma::vec rhs_AlphaZ =
      (1.0 / sigma_Z_sq) * V.t() * (z - H * Pi_Z.t()) +
      (1.0 / sigma_Y_sq) * VWz.t() * (y - mu
        - X_Dash * Beta_X
        - V * Alpha_Y.t()
        - (H * Pi_Z.t()) % WZ_diag);

    arma::rowvec m_Alpha_Z = (Sigma_Alpha_Z * rhs_AlphaZ).t();
    Alpha_Z = arma::mvnrnd(m_Alpha_Z.t(), Sigma_Alpha_Z, 1).t();

    // ---- Step 6: sigma_Z^2 | rest -------------------------------------------
    arma::vec z_res = z - (H * Pi_Z.t() + V * Alpha_Z.t());
    double penZ = arma::as_scalar(Pi_Z * H_Star * Pi_Z.t()) / gH
      + arma::as_scalar(Alpha_Z * V_Star * Alpha_Z.t()) / gV;
    double sZ_shape = a_sigma + 0.5 * (n + pH + pV);
    double sZ_scale = b_sigma + 0.5 * (arma::dot(z_res, z_res) + penZ);
    sigma_Z_sq = rinvgamma_1(sZ_shape, sZ_scale);

    // =========================================================================
    // BLOCK 3 -- Outcome second stage: mu, Alpha_Y, Beta_Z, sigma_Y^2,
    //            (Beta_X, Beta_XZ), gamma, rho
    //
    // Recompute IV-projected values using freshly updated first-stage draws
    // before entering any second-stage update.
    // =========================================================================
    X_Dash              = G * Pi_X.t() + V * Alpha_X.t();  // X*      (n x 1)
    Z_Dash              = H * Pi_Z.t() + V * Alpha_Z.t();  // Z*      (n x 1)
    arma::vec XZ_Dash   = X_Dash % Z_Dash;                  // X* o Z* (n x 1)

    // Design matrix for the causal block [X*, X*Z*]
    arma::mat XBeta_Dash(n, 2);
    XBeta_Dash.col(0) = X_Dash;
    XBeta_Dash.col(1) = XZ_Dash;

    arma::mat XBeta_Star = XBeta_Dash.t() * XBeta_Dash;   // (2 x 2)

    // Predicted Y from causal block -- used as offset in downstream steps
    arma::vec XB_Beta = X_Dash * Beta_X + XZ_Dash * Beta_XZ;                     // (n x 1)

    // ---- Step 7: mu | rest --------------------------------------------------
    // Prior:     mu ~ N(0, c_mu * sigma_Y^2)
    // Posterior: mu | rest ~ N(m_mu, v_mu)
    //   v_mu = sigma_Y^2 / (n + 1/c_mu)
    //   m_mu = sum(r_mu)  / (n + 1/c_mu),  r_mu = y - XB_Beta - Z*Beta_Z - V*Alpha_Y
    arma::vec r_mu = y - XB_Beta - Z_Dash * Beta_Z - V * Alpha_Y.t();
    double v_mu = sigma_Y_sq / (static_cast<double>(n) + 1.0 / c_mu);
    double m_mu = arma::sum(r_mu) / (static_cast<double>(n) + 1.0 / c_mu);
    mu = R::rnorm(m_mu, std::sqrt(v_mu));

    // ---- Step 8: Alpha_Y | rest ---------------------------------------------
    // Zellner g-prior posterior: shrinkage factor cV = gV / (1 + gV).
    double cV = gV / (1.0 + gV);
    arma::vec    y_resY    = y - mu - XB_Beta - Z_Dash * Beta_Z;
    arma::rowvec m_Alpha_Y = (cV * V_Tilde * y_resY).t();
    arma::mat    S_Alpha_Y = (cV * sigma_Y_sq) * V_Star_Inv;
    Alpha_Y = arma::mvnrnd(m_Alpha_Y.t(), S_Alpha_Y, 1).t();

    // ---- Step 9: Beta_Z | rest ----------------------------------------------
    // Scalar conjugate normal update. Prior Beta_Z ~ N(0, gZ sigma_Y^2 / d_Z)
    // with d_Z fixed, so the posterior precision is (Z*'Z* + d_Z/gZ)/sigma_Y^2
    // and the mean is Z*'r / (Z*'Z* + d_Z/gZ).
    double ZtZ       = arma::dot(Z_Dash, Z_Dash);
    arma::vec y_resZ = y - mu - XB_Beta - V * Alpha_Y.t();
    double aZ        = ZtZ + d_Z / gZ;
    double m_BetaZ   = arma::dot(Z_Dash, y_resZ) / aZ;
    double v_BetaZ   = sigma_Y_sq / aZ;
    Beta_Z = R::rnorm(m_BetaZ, std::sqrt(v_BetaZ));

    // ---- Step 10: sigma_Y^2 | rest ------------------------------------------
    // Conjugate InvGamma. The scale includes g-prior quadratic penalty terms for
    // the intercept (mu), the causal block (Beta_X, Beta_XZ), the receptor effect
    // (Beta_Z), and the covariate block (Alpha_Y).
    // Degrees contributed: n (data) + pV (Alpha_Y) + 4 (mu, Beta_X, Beta_XZ, Beta_Z).
    arma::vec resY    = y - mu - XB_Beta - Z_Dash * Beta_Z - V * Alpha_Y.t();
    double s_gamma    = (gamma == 1) ? 1.0 : nu1;  // slab or spike variance scale

    // Prior precision shape for the causal block (the fixed D0) and for
    // Beta_Z (the fixed d_Z). Used here for the sigma_Y^2 penalty and again in
    // Steps 11 and 12.
    const arma::mat& B_shape = D0;
    const double     z_shape = d_Z;

    double pen_mu = (mu * mu) / c_mu;
    double pen_B  = arma::as_scalar(Beta_row * B_shape * Beta_row.t()) / (gBeta * s_gamma);
    double pen_BZ = (Beta_Z * Beta_Z) * z_shape / gZ;
    double pen_AY = arma::as_scalar(Alpha_Y * V_Star * Alpha_Y.t()) / gV;

    double sY_shape = a_sigma + 0.5 * (n + pV + 4.0);
    double sY_scale = b_sigma + 0.5 * (arma::dot(resY, resY)
      + pen_mu + pen_B + pen_BZ + pen_AY);
    sigma_Y_sq = rinvgamma_1(sY_shape, sY_scale);

    // ---- Step 11: (Beta_X, Beta_XZ) | rest  [spike-and-slab] ---------------
    // Prior N(0, gBeta s_gamma sigma_Y^2 D0^{-1}), slab (s_gamma = 1) or spike
    // (s_gamma = nu1). Normal prior, normal likelihood, so the posterior is
    //   N( A^{-1} X_beta' r,  sigma_Y^2 A^{-1} ),  A = X_beta'X_beta + D0/(gBeta s_gamma).
    // One 2 x 2 solve.
    arma::vec y_res_B = y - mu - Z_Dash * Beta_Z - V * Alpha_Y.t();
    arma::vec rhsB    = XBeta_Dash.t() * y_res_B;  // sufficient statistic (2 x 1)

    {
      arma::vec mB(2);
      arma::mat SB(2, 2);
      arma::mat A_inv = safe_inv_sympd(
        XBeta_Star + D0 / (gBeta * s_gamma) + ridge * arma::eye(2, 2),
        "causal-block precision");
      mB = A_inv * rhsB;
      SB = sigma_Y_sq * A_inv;
      arma::vec draw_b = arma::mvnrnd(mB, SB, 1);
      Beta_X  = draw_b(0);
      Beta_XZ = draw_b(1);
    }
    Beta_row(0) = Beta_X;
    Beta_row(1) = Beta_XZ;

    // ---- Step 12: gamma | rest  [Bernoulli spike-and-slab] ------------------
    // Posterior probability of gamma = 1 from the ratio of the two prior
    // densities of the current Beta draw (slab vs. spike) and the prior on
    // rho. The -log(nu1) in logB is the two-dimensional normalising-constant
    // difference; the quadratic form uses the same shape as the prior.
    double quad = arma::as_scalar(Beta_row * B_shape * Beta_row.t());
    double logA = -0.5 * quad / (sigma_Y_sq * gBeta)
      + std::log(rho);
    double logB = -0.5 * quad / (sigma_Y_sq * gBeta * nu1)
      + std::log(1.0 - rho) - std::log(nu1);
    double lse  = logsumexp2(logA, logB);
    double p    = std::exp(logA - lse);

    // Guard against numerical under-/overflow
    if (!R_finite(p))      p = 0.5;
    if (p < 1e-12)         p = 1e-12;
    if (p > 1.0 - 1e-12)  p = 1.0 - 1e-12;
    gamma = R::rbinom(1.0, p);

    // ---- Step 13: rho | rest ------------------------------------------------
    // Conjugate Beta update:  rho | gamma ~ Beta(a_rho + gamma, b_rho + 1 - gamma)
    rho = R::rbeta(a_rho + gamma, b_rho + 1 - gamma);

    // =========================================================================
    // Accumulate posterior-mean estimates (post-burn-in, thinned)
    // =========================================================================
    if (it > burn_in && ((it - burn_in) % thin == 0)) {
      ++save_idx;
      Pi_X_sum    += Pi_X;
      Alpha_X_sum += Alpha_X;
      Pi_Z_sum    += Pi_Z;
      Alpha_Z_sum += Alpha_Z;
      Beta_X_sum  += Beta_X;
      Beta_XZ_sum += Beta_XZ;
      Beta_Z_sum  += Beta_Z;
      Alpha_Y_sum += Alpha_Y;
      mu_sum      += mu;
      sX_sum      += sigma_X_sq;
      sZ_sum      += sigma_Z_sq;
      sY_sum      += sigma_Y_sq;
      gamma_sum   += static_cast<double>(gamma);
      rho_sum     += rho;

      // Store the scalar draws (guarded against any off-by-one in n_keep)
      if (save_idx <= n_keep) {
        const int si = save_idx - 1;
        Beta_X_draws(si)  = Beta_X;
        Beta_XZ_draws(si) = Beta_XZ;
        Beta_Z_draws(si)  = Beta_Z;
        gamma_draws(si)   = static_cast<double>(gamma);
        mu_draws(si)      = mu;

        // ---- Joint log-likelihood at the current draw ------------------------
        // The three equations are conditionally independent given the
        // parameters, so the joint log-likelihood is the sum of three Gaussian
        // log-densities:
        //
        //   X | Pi_X, Alpha_X, sigma_X^2 ~ N(G Pi_X' + V Alpha_X', sigma_X^2 I)
        //   Z | Pi_Z, Alpha_Z, sigma_Z^2 ~ N(H Pi_Z' + V Alpha_Z', sigma_Z^2 I)
        //   Y | ...                      ~ N(mu + X* Beta_X + X*Z* Beta_XZ
        //                                      + Z* Beta_Z + V Alpha_Y',
        //                                    sigma_Y^2 I)
        //
        // x_res and z_res were formed in Steps 3 and 6 from Pi_X, Alpha_X,
        // Pi_Z and Alpha_Z, none of which are revisited later in the sweep, so
        // they remain the residuals of the current draw. X_Dash, Z_Dash and
        // XZ_Dash (Step 7) are current for the same reason. Only the causal
        // block was redrawn afterwards, in Step 11, so the Y-equation fitted
        // value is recomputed here from the final Beta_X and Beta_XZ.
        //
        // This is the likelihood alone: the g-prior quadratic penalties that
        // enter the variance updates are deliberately excluded, so the quantity
        // is a measure of fit and not an unnormalised posterior.
        const arma::vec XB_Beta_cur = X_Dash * Beta_X + XZ_Dash * Beta_XZ;
        const arma::vec resY_cur    = y - mu - XB_Beta_cur
                                        - Z_Dash * Beta_Z - V * Alpha_Y.t();

        const double nd = static_cast<double>(n);
        loglik_draws(si) =
          -0.5 * nd * (LOG_2PI + std::log(sigma_X_sq))
          - 0.5 * arma::dot(x_res, x_res) / sigma_X_sq
          - 0.5 * nd * (LOG_2PI + std::log(sigma_Z_sq))
          - 0.5 * arma::dot(z_res, z_res) / sigma_Z_sq
          - 0.5 * nd * (LOG_2PI + std::log(sigma_Y_sq))
          - 0.5 * arma::dot(resY_cur, resY_cur) / sigma_Y_sq;
      }
    }

    // Progress reporting is off by default: a full analysis runs several
    // chains on each of many triplets.
    if (verbose && (it % 1000 == 0))
      Rcpp::Rcout << "Iteration " << it << " / " << n_iter << " completed.\n";

  }  // end MCMC loop

  // --------------------------------------------------------------------------
  // Return posterior means
  // --------------------------------------------------------------------------
  const double denom = static_cast<double>(save_idx > 0 ? save_idx : 1);

  // Rcpp's List::create accepts at most twenty arguments, and the return value
  // carries twenty-two elements. The list is therefore allocated at its final
  // length and filled position by position, with the names assigned together at
  // the end so that each name stays adjacent to the value it labels.
  const int n_out = 22;
  List out(n_out);
  CharacterVector out_names(n_out);
  int k = 0;

  // ---- Posterior means -----------------------------------------------------
  out[k] = Pi_X_sum    / denom; out_names[k] = "Pi_X_mean";       ++k;  // eQTL effects on ligand (sender)
  out[k] = Alpha_X_sum / denom; out_names[k] = "Alpha_X_mean";    ++k;  // covariate effects on ligand
  out[k] = Pi_Z_sum    / denom; out_names[k] = "Pi_Z_mean";       ++k;  // eQTL effects on receptor (receiver)
  out[k] = Alpha_Z_sum / denom; out_names[k] = "Alpha_Z_mean";    ++k;  // covariate effects on receptor
  out[k] = Beta_X_sum  / denom; out_names[k] = "Beta_X_mean";     ++k;  // causal effect: ligand -> pathway
  out[k] = Beta_XZ_sum / denom; out_names[k] = "Beta_XZ_mean";    ++k;  // interaction effect (receptor-modulated)
  out[k] = Beta_Z_sum  / denom; out_names[k] = "Beta_Z_mean";     ++k;  // direct receptor effect on pathway
  out[k] = Alpha_Y_sum / denom; out_names[k] = "Alpha_Y_mean";    ++k;  // covariate effects on pathway
  out[k] = mu_sum      / denom; out_names[k] = "mu_mean";         ++k;  // intercept
  out[k] = sX_sum      / denom; out_names[k] = "sigma_X_sq_mean"; ++k;  // sender variance
  out[k] = sZ_sum      / denom; out_names[k] = "sigma_Z_sq_mean"; ++k;  // receiver variance
  out[k] = sY_sum      / denom; out_names[k] = "sigma_Y_sq_mean"; ++k;  // outcome variance
  out[k] = gamma_sum   / denom; out_names[k] = "gamma_mean";      ++k;  // posterior inclusion probability (PIP)
  out[k] = rho_sum     / denom; out_names[k] = "rho_mean";        ++k;  // posterior mean inclusion prior

  // ---- Retained draws (length n_keep) --------------------------------------
  // Needed for credible intervals on functionals such as
  // tau = -Beta_X/Beta_XZ, for Monte Carlo standard errors and effective
  // sample size on the PIP, and for trace plots.
  out[k] = save_idx;                        out_names[k] = "n_keep";        ++k;
  out[k] = Beta_X_draws.head(save_idx);     out_names[k] = "Beta_X_draws";  ++k;
  out[k] = Beta_XZ_draws.head(save_idx);    out_names[k] = "Beta_XZ_draws"; ++k;
  out[k] = Beta_Z_draws.head(save_idx);     out_names[k] = "Beta_Z_draws";  ++k;
  out[k] = gamma_draws.head(save_idx);      out_names[k] = "gamma_draws";   ++k;
  out[k] = mu_draws.head(save_idx);         out_names[k] = "mu_draws";      ++k;
  out[k] = loglik_draws.head(save_idx);     out_names[k] = "loglik_draws";  ++k;

  // ---- Prior scale actually used ------------------------------------------
  // Recorded so that any fit can be checked against the scale it was run
  // under.
  NumericVector prior_scale = NumericVector::create(
    Named("d_X")  = d_X,
    Named("d_XZ") = d_XZ,
    Named("d_Z")  = d_Z);
  out[k] = prior_scale;                     out_names[k] = "prior_scale";   ++k;

  out.attr("names") = out_names;
  return out;
}
