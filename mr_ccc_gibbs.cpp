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
//    Y = mu  +  gamma * (X* * Beta_X  +  X*Z* * Beta_XZ)
//           +  Z* * Beta_Z  +  V * Alpha_Y  +  eps_Y,
//                                               eps_Y ~ N(0, sigma_Y^2 * I)
//
//  where X* = G*Pi_X + V*Alpha_X  and  Z* = H*Pi_Z + V*Alpha_Z
//  are the IV-projected (genetically predicted) expression levels.
//
//  Spike-and-slab prior on the causal block (Beta_X, Beta_XZ):
//    gamma = 1 (slab):  (Beta_X, Beta_XZ) ~ N(0, gBeta * sigma_Y^2 * ([X* X*Z*]'[X* X*Z*])^{-1})
//    gamma = 0 (spike): same form but variance scaled down by nu1  (nu1 << 1)
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
//  Named R list of posterior means for all model parameters.
// =============================================================================

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


// =============================================================================
// Main exported function
// =============================================================================

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
    int    thin     = 10,     // thinning interval
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
    double ridge    = 1e-8    // diagonal ridge for stable matrix inversion
) {

  // --------------------------------------------------------------------------
  // Dimensions
  // --------------------------------------------------------------------------
  const int n  = X.n_rows;
  const int pG = G.n_cols;
  const int pH = H.n_cols;
  const int pV = V.n_cols;

  // --------------------------------------------------------------------------
  // Precompute Gram matrices and regularized inverses.
  // These are fixed throughout the MCMC and reused at every iteration.
  // --------------------------------------------------------------------------
  const arma::mat G_Star = G.t() * G;  // (pG x pG)
  const arma::mat V_Star = V.t() * V;  // (pV x pV)
  const arma::mat H_Star = H.t() * H;  // (pH x pH)

  const arma::mat G_Star_Inv = arma::inv_sympd(G_Star + ridge * arma::eye(pG, pG));
  const arma::mat V_Star_Inv = arma::inv_sympd(V_Star + ridge * arma::eye(pV, pV));
  const arma::mat H_Star_Inv = arma::inv_sympd(H_Star + ridge * arma::eye(pH, pH));

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
  // Column vectors extracted once (avoids repeated .col(0) inside the loop)
  // --------------------------------------------------------------------------
  const arma::vec x = X.col(0);
  const arma::vec z = Z.col(0);
  const arma::vec y = Y.col(0);

  // ==========================================================================
  // MCMC (Gibbs) loop
  // ==========================================================================
  int save_idx = 0;

  for (int it = 1; it <= n_iter; ++it) {

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

    arma::mat Sigma_Pi_X = arma::inv_sympd(
      (1.0 / sigma_X_sq) * (1.0 + 1.0 / gG) * G_Star +
      (1.0 / sigma_Y_sq) * (GW.t() * GW) +
      ridge * arma::eye(pG, pG));

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

    arma::mat Sigma_Alpha_X = arma::inv_sympd(
      (1.0 / sigma_X_sq) * (1.0 + 1.0 / gV) * V_Star +
      (1.0 / sigma_Y_sq) * (VW.t() * VW) +
      ridge * arma::eye(pV, pV));

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

    arma::mat Sigma_Pi_Z = arma::inv_sympd(
      (1.0 / sigma_Z_sq) * (1.0 + 1.0 / gH) * H_Star +
      (1.0 / sigma_Y_sq) * (HW.t() * HW) +
      ridge * arma::eye(pH, pH));

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

    arma::mat Sigma_Alpha_Z = arma::inv_sympd(
      (1.0 / sigma_Z_sq) * (1.0 + 1.0 / gV) * V_Star +
      (1.0 / sigma_Y_sq) * (VWz.t() * VWz) +
      ridge * arma::eye(pV, pV));

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

    arma::mat XBeta_Star     = XBeta_Dash.t() * XBeta_Dash;                      // (2 x 2)
    arma::mat XBeta_Star_Inv = arma::inv_sympd(XBeta_Star + ridge * arma::eye(2, 2));

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
    // Scalar conjugate normal update under the g-prior for Beta_Z.
    double cZ        = gZ / (1.0 + gZ);
    double ZtZ       = arma::dot(Z_Dash, Z_Dash);
    arma::vec y_resZ = y - mu - XB_Beta - V * Alpha_Y.t();
    double m_BetaZ   = cZ * arma::dot(Z_Dash, y_resZ) / (ZtZ + ridge);
    double v_BetaZ   = cZ * sigma_Y_sq / (ZtZ + ridge);
    Beta_Z = R::rnorm(m_BetaZ, std::sqrt(v_BetaZ));

    // ---- Step 10: sigma_Y^2 | rest ------------------------------------------
    // Conjugate InvGamma. The scale includes g-prior quadratic penalty terms for
    // the intercept (mu), the causal block (Beta_X, Beta_XZ), the receptor effect
    // (Beta_Z), and the covariate block (Alpha_Y).
    // Degrees contributed: n (data) + pV (Alpha_Y) + 4 (mu, Beta_X, Beta_XZ, Beta_Z).
    arma::vec resY    = y - mu - XB_Beta - Z_Dash * Beta_Z - V * Alpha_Y.t();
    double s_gamma    = (gamma == 1) ? 1.0 : nu1;  // slab or spike variance scale
    arma::mat S_B_inv = (1.0 / gBeta) * XBeta_Star;

    double pen_mu = (mu * mu) / c_mu;
    double pen_B  = arma::as_scalar(Beta_row * S_B_inv * Beta_row.t()) / s_gamma;
    double pen_BZ = (Beta_Z * Beta_Z) * ZtZ / gZ;
    double pen_AY = arma::as_scalar(Alpha_Y * V_Star * Alpha_Y.t()) / gV;

    double sY_shape = a_sigma + 0.5 * (n + pV + 4.0);
    double sY_scale = b_sigma + 0.5 * (arma::dot(resY, resY)
      + pen_mu + pen_B + pen_BZ + pen_AY);
    sigma_Y_sq = rinvgamma_1(sY_shape, sY_scale);

    // ---- Step 11: (Beta_X, Beta_XZ) | rest  [spike-and-slab] ---------------
    // When gamma = 1 (slab):  g-prior with effective scale gBeta.
    // When gamma = 0 (spike): g-prior with scale nu1 * gBeta  (near-zero prior).
    // Both cases share the same MVN form; only the effective g-scale differs.
    arma::vec y_res_B = y - mu - Z_Dash * Beta_Z - V * Alpha_Y.t();
    arma::vec rhsB    = XBeta_Dash.t() * y_res_B;  // sufficient statistic (2 x 1)

    {
      double g_eff     = (gamma == 1) ? gBeta : (nu1 * gBeta);
      double cB        = g_eff / (1.0 + g_eff);
      arma::vec  mB    = cB * (XBeta_Star_Inv * rhsB);
      arma::mat  SB    = cB * sigma_Y_sq * XBeta_Star_Inv;
      arma::vec draw_b = arma::mvnrnd(mB, SB, 1);
      Beta_X  = draw_b(0);
      Beta_XZ = draw_b(1);
    }
    Beta_row(0) = Beta_X;
    Beta_row(1) = Beta_XZ;

    // ---- Step 12: gamma | rest  [Bernoulli spike-and-slab] ------------------
    // Posterior probability of gamma = 1 via the marginal likelihood ratio
    // (slab vs. spike) combined with the Beta prior on rho.
    double quad = arma::as_scalar(Beta_row * XBeta_Star * Beta_row.t());
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
    }

    if (it % 1000 == 0)
      Rcpp::Rcout << "Iteration " << it << " / " << n_iter << " completed.\n";

  }  // end MCMC loop

  // --------------------------------------------------------------------------
  // Return posterior means
  // --------------------------------------------------------------------------
  const double denom = static_cast<double>(save_idx > 0 ? save_idx : 1);

  return List::create(
    _["Pi_X_mean"]       = Pi_X_sum    / denom,  // eQTL effects on ligand (sender)
    _["Alpha_X_mean"]    = Alpha_X_sum / denom,  // covariate effects on ligand
    _["Pi_Z_mean"]       = Pi_Z_sum    / denom,  // eQTL effects on receptor (receiver)
    _["Alpha_Z_mean"]    = Alpha_Z_sum / denom,  // covariate effects on receptor
    _["Beta_X_mean"]     = Beta_X_sum  / denom,  // causal effect: ligand -> pathway
    _["Beta_XZ_mean"]    = Beta_XZ_sum / denom,  // interaction effect (receptor-modulated)
    _["Beta_Z_mean"]     = Beta_Z_sum  / denom,  // direct receptor effect on pathway
    _["Alpha_Y_mean"]    = Alpha_Y_sum / denom,  // covariate effects on pathway
    _["mu_mean"]         = mu_sum      / denom,  // intercept
    _["sigma_X_sq_mean"] = sX_sum      / denom,  // sender variance
    _["sigma_Z_sq_mean"] = sZ_sum      / denom,  // receiver variance
    _["sigma_Y_sq_mean"] = sY_sum      / denom,  // outcome variance
    _["gamma_mean"]      = gamma_sum   / denom,  // posterior inclusion probability (PIP)
    _["rho_mean"]        = rho_sum     / denom   // posterior mean inclusion prior
  );
}
