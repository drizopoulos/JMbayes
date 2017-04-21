#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

uvec seqC(int x, int y, int by) {
    int nRatio = (y - x) / by;
    uvec anOut(nRatio + 1);
    int n = 0;
    for (int i = x; i <= y; i = i + by) {
        anOut[n] = i;
        n += 1;
    }
    return anOut;
}

arma::vec rowsum(const vec& x_, const uvec& group) {
    vec x = cumsum(x_);
    vec out = x.elem(group);
    out.insert_rows(0, 1);
    unsigned int n = out.n_elem;
    out = out.rows(1, n - 1) - out.rows(0, n - 2);
    return(out);
}

arma::vec Vpnorm(const vec& x) {
    int n = x.size();
    vec res(n);
    for(int i = 0; i < n; ++i) {
        res[i] = ::Rf_pnorm5(x[i], 0.0, 1.0, 1, 0);
    }
    return(res);
}

arma::field<arma::vec> List2Field_vec (const List & Vecs) {
    int n_list = Vecs.size();
    arma::field<arma::vec> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<arma::vec>(Vecs[i]);
    }
    return(res);
}

arma::field<arma::uvec> List2Field_uvec (const List & uVecs, bool substract1 = true) {
    int n_list = uVecs.size();
    arma::field<arma::uvec> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        if (substract1) {
            res.at(i) = as<arma::uvec>(uVecs[i]) - 1;
        } else {
            res.at(i) = as<arma::uvec>(uVecs[i]);
        }
    }
    return(res);
}

arma::field<arma::mat> List2Field_mat (const List & Mats) {
    int n_list = Mats.size();
    arma::field<arma::mat> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<arma::mat>(Mats[i]);
    }
    return(res);
}

arma::field<arma::vec> lin_predF(const field<vec>& Xbetas, const field<mat>& Z,
                                 const mat& b, const field<uvec>& RE_inds,
                                 const field<uvec>& id) {
    signed int n_field = Xbetas.size();
    field<vec> out(n_field);
    for (int i = 0; i < n_field; ++i) {
        mat bb = b.cols(RE_inds.at(i));
        vec Zb = sum(Z.at(i) % bb.rows(id.at(i)), 1);
        out.at(i) = Xbetas.at(i) + Zb;
    }
    return(out);
}

arma::mat lin_pred_matF(const field<vec>& Xbetas, const field<mat>& Z,
                        const mat& b, const field<mat>& U,
                        const field<uvec>& RE_inds, const field<uvec>& id,
                        const field<uvec>& col_inds, const uvec& row_inds,
                        const int& nrows, const int& ncols) {
    int n_field = Xbetas.size();
    mat out = mat(nrows, ncols, fill::zeros);
    for (int i = 0; i < n_field; ++i) {
        mat bb = b.cols(RE_inds.at(i));
        vec Zb = sum(Z.at(i) % bb.rows(id.at(i)), 1);
        out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % (Xbetas.at(i) + Zb);
    }
    return(out);
}

mat mvrnorm(const int& n, const mat& var_cov, bool cholsigma = false) {
    int ncols = var_cov.n_cols;
    mat Y = randn(n, ncols);
    if (cholsigma == false) {
        mat H = chol(var_cov);
        Y = trans(Y * H);
    } else {
        Y = trans(Y * var_cov);
    }
    return(Y);
}

cube mvrnorm_cube(const int& n, const cube& var_cov, vec sigmas, bool cholsigma = true) {
    unsigned int k = var_cov.n_rows;
    unsigned int N = var_cov.n_slices;
    cube out = cube(n, k, N, fill::zeros);
    for (uword i = 0; i < N; ++i) {
        if (cholsigma == true) {
            out.slice(i) = sqrt(sigmas[i]) * (randn(n, k) * var_cov.slice(i));
        } else {
            mat H = chol(var_cov.slice(i));
            out.slice(i) = sqrt(sigmas[i]) * (randn(n, k) * H);
        }
    }
    return(out);
}

arma::mat extract_b(const cube& b, unsigned int sample) {
    unsigned int n = b.n_slices;
    unsigned int k = b.n_cols;
    mat out = mat(n, k, fill::zeros);
    for (uword i = 0; i < n; ++i) {
        out.row(i) = b.slice(i).row(sample);
    }
    return(out);
}

arma::vec log_p_event_RC (const arma::vec& logh, const arma::vec& H, 
                          const arma::vec& event) {
    arma::vec res = event % logh - H;
    return(res);
}

arma::vec log_p_event_IC (const arma::vec& logh, const arma::vec& H, const arma::vec& HU,
                          const LogicalVector& Levent1, const LogicalVector& Levent01,
                          const LogicalVector& Levent2, const LogicalVector& Levent3) {
    int n = logh.n_rows;
    arma::vec res(n, fill::zeros);
    for (int i = 0; i < n; ++i) {
        if (Levent1(i)) res.at(i) += logh.at(i);
        if (Levent01(i)) res.at(i) -= H.at(i);
        if (Levent2(i)) res.at(i) += log(1 - exp(-H.at(i)));
        if (Levent3(i)) res.at(i) += log(exp(-HU.at(i)) - exp(HU.at(i)));
    }
    return(res);
}

arma::vec log_longF(const field<vec>& y, const field<vec>& eta,
                    const CharacterVector& fams, const CharacterVector& links,
                    const List& sigmas, const field<uvec>& id, const int& n) {
    int n_outcomes = y.size();
    vec out(n, fill::zeros);
    for (int i = 0; i < n_outcomes; ++i) {
        uvec id_i = id.at(i);
        vec y_i = y.at(i);
        vec eta_i = eta.at(i);
        if (fams[i] == "gaussian") {
            double sigma_i = as<double>(sigmas[i]);
            vec log_dens = - 0.5 * pow((y_i - eta_i) / sigma_i, 2);
            out += rowsum(log_dens, id_i);
        } else if (fams[i] == "binomial") {
            if (links[i] == "logit") {
                vec pr = exp(eta_i) / (1 + exp(eta_i));
                vec log_dens = y_i % log(pr) + (1 - y_i) % log(1 - pr);
                out += rowsum(log_dens, id_i);
            } else if (links[i] == "probit") {
                vec pr = Vpnorm(eta_i);
                vec log_dens = y_i % log(pr) + (1 - y_i) % log(1 - pr);
                out += rowsum(log_dens, id_i);
            } else if (links[i] == "cloglog") {
                vec pr = - exp(- exp(eta_i)) + 1;
                vec log_dens = y_i % log(pr) + (1 - y_i) % log(1 - pr);
                out += rowsum(log_dens, id_i);
            }
        } else if (fams[i] == "poisson") {
            vec mu = exp(eta_i);
            vec log_dens = y_i % log(mu) - mu;
            out += rowsum(log_dens, id_i);
        }
    }
    return(out);
}

arma::mat rwishart(int df, const arma::mat& S) {
    int m = S.n_rows;
    mat Z(m, m);
    for (int i = 0; i < m; ++i) {
        Z(i, i) = sqrt(R::rchisq(df - i));
    }
    for (int j = 0; j < m; ++j){
        for (int i = j + 1; i < m; ++i){
            Z(i, j) = R::rnorm(0, 1);
        }
    }
    arma::mat C = arma::trimatl(Z).t() * arma::chol(S);
    return C.t() * C;
}

arma::vec log_postREF(const mat& b_mat, const vec& Bs_gammas, const vec& gammas,
                      const vec& alphas, const field<vec>& y, const field<vec>& Xbetas,
                      const field<mat>& Z,
                      const field<uvec>& RE_inds, const field<uvec>& RE_inds2,
                      const field<uvec>& idL, const field<uvec>& idL2,
                      const CharacterVector& fams, const CharacterVector& links,
                      const List& sigmas, const mat& invD,
                      const int& n, const int& ns, const int& n_alphas,
                      const vec& event, const mat& W1, const mat& W1s, const mat& W2,
                      const mat& W2s, const field<vec>& XXbetas,
                      const field<vec>& XXsbetas, const field<mat>& ZZ,
                      const field<mat>& ZZs, const field<mat>& U,
                      const field<mat>& Us, const vec& Pw, const uvec& idGK,
                      const field<uvec>& idT, const field<uvec>& idTs,
                      const field<uvec>& col_inds,
                      const uvec& row_inds_U, const uvec& row_inds_Us) {
    field<vec> eta = lin_predF(Xbetas, Z, b_mat, RE_inds, idL);
    vec log_pyb = log_longF(y, eta, fams, links, sigmas, idL2, n);
    vec log_pb = - 0.5 * sum((b_mat * invD) % b_mat, 1);
    mat Wlong = lin_pred_matF(XXbetas, ZZ, b_mat, U, RE_inds2, idT,
                              col_inds, row_inds_U, n, n_alphas);
    mat Wlongs = lin_pred_matF(XXsbetas, ZZs, b_mat, Us, RE_inds2, idTs,
                               col_inds, row_inds_Us, ns, n_alphas);
    vec log_h = W1 * Bs_gammas + W2 * gammas + Wlong * alphas;
    vec H = rowsum(Pw % exp(W1s * Bs_gammas + W2s * gammas + Wlongs * alphas), idGK);
    vec log_ptb = (event % log_h) - H;
    vec out = log_pyb + log_ptb + log_pb;
    return(out);
}

arma::vec log_postRE_nogammasF(const mat& b_mat, const vec& Bs_gammas, const vec& alphas,
                               const field<vec>& y, const field<vec>& Xbetas,
                               const field<mat>& Z,
                               const field<uvec>& RE_inds, const field<uvec>& RE_inds2,
                               const field<uvec>& idL, const field<uvec>& idL2,
                               const CharacterVector& fams, const CharacterVector& links,
                               const List& sigmas, const mat& invD,
                               const int& n, const int& ns, const int& n_alphas,
                               const vec& event, const mat& W1,
                               const mat& W1s, const field<vec>& XXbetas,
                               const field<vec>& XXsbetas, const field<mat>& ZZ,
                               const field<mat>& ZZs, const field<mat>& U,
                               const field<mat>& Us, const vec& Pw, const uvec& idGK,
                               const field<uvec>& idT, const field<uvec>& idTs,
                               const field<uvec>& col_inds,
                               const uvec& row_inds_U, const uvec& row_inds_Us) {
    field<vec> eta = lin_predF(Xbetas, Z, b_mat, RE_inds, idL);
    vec log_pyb = log_longF(y, eta, fams, links, sigmas, idL2, n);
    vec log_pb = - 0.5 * sum((b_mat * invD) % b_mat, 1);
    mat Wlong = lin_pred_matF(XXbetas, ZZ, b_mat, U, RE_inds2, idT,
                              col_inds, row_inds_U, n, n_alphas);
    mat Wlongs = lin_pred_matF(XXsbetas, ZZs, b_mat, Us, RE_inds2, idTs,
                               col_inds, row_inds_Us, ns, n_alphas);
    vec log_h = W1 * Bs_gammas + Wlong * alphas;
    vec H = rowsum(Pw % exp(W1s * Bs_gammas + Wlongs * alphas), idGK);
    vec log_ptb = (event % log_h) - H;
    vec out = log_pyb + log_ptb + log_pb;
    return(out);
}

double log_weightsREF(const mat& b_mat, const field<vec>& y, const field<vec>& Xbetas,
                      const field<mat>& Z, const field<uvec>& RE_inds,
                      const field<uvec>& idL, const field<uvec>& idL2,
                      const CharacterVector& fams, const CharacterVector& links,
                      const List& sigmas, const mat& invD, const int& n) {
    field<vec> eta = lin_predF(Xbetas, Z, b_mat, RE_inds, idL);
    vec log_pyb = log_longF(y, eta, fams, links, sigmas, idL2, n);
    vec log_pb = - 0.5 * sum((b_mat * invD) % b_mat, 1);
    double out = - sum(log_pyb + log_pb);
    return(out);
}

double logPrior(const vec& x, const vec& mean, const mat& Tau, double tau = 1) {
    vec z = (x - mean);
    return(- 0.5 * tau * as_scalar(z.t() * Tau * z));
}

double logPosterior (const vec& event,
                     const mat& W1, const mat& W1s, const vec& Bs_gammas,
                     const mat& W2, const mat& W2s, const vec& gammas,
                     const mat& Wlong, const mat& Wlongs, const vec& alphas,
                     const vec& Pw,
                     const vec& mean_Bs_gammas, const mat& Tau_Bs_gammas, double tau_Bs_gammas,
                     const vec& mean_gammas, const mat& Tau_gammas, double tau_gammas,
                     const vec& mean_alphas, const mat& Tau_alphas, double tau_alphas) {
    // log-likelihood
    vec log_h = W1 * Bs_gammas + W2 * gammas + Wlong * alphas;
    vec H = Pw % exp(W1s * Bs_gammas + W2s * gammas + Wlongs * alphas);
    //
    double out = sum(event % log_h) - sum(H) +
        logPrior(Bs_gammas, mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas) +
        logPrior(gammas, mean_gammas, Tau_gammas, tau_gammas) +
        logPrior(alphas, mean_alphas, Tau_alphas, tau_alphas);
    return(out);
}

double logPosterior_nogammas (const vec& event,
                     const mat& W1, const mat& W1s, const vec& Bs_gammas,
                     const mat& Wlong, const mat& Wlongs, const vec& alphas,
                     const vec& Pw,
                     const vec& mean_Bs_gammas, const mat& Tau_Bs_gammas, double tau_Bs_gammas,
                     const vec& mean_alphas, const mat& Tau_alphas, double tau_alphas) {
    // log-likelihood
    vec log_h = W1 * Bs_gammas + Wlong * alphas;
    vec H = Pw % exp(W1s * Bs_gammas + Wlongs * alphas);
    //
    double out = sum(event % log_h) - sum(H) +
        logPrior(Bs_gammas, mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas) +
        logPrior(alphas, mean_alphas, Tau_alphas, tau_alphas);
    return(out);
}

double bounds_sigma (double sigma, double eps1, double eps3) {
    if (sigma < eps1) {
        sigma = eps1;
    }
    if (sigma > eps3) {
        sigma = eps3;
    }
    return(sigma);
}

mat bounds_Cov (mat Sigma, double eps2, double eps3) {
    int n = Sigma.n_rows;
    double D = det(Sigma);
    if (D > eps3) {
        Sigma = (eps3 / D) * Sigma;
    }
    return(Sigma + eps2 * diagmat(vec(n, fill::ones)));
}

// [[Rcpp::export]]
List lap_rwm_C (List initials, List Data, List priors, List scales, List Covs,
                List control, bool interval_cens) {
    // Data
    List y = as<List>(Data["y"]);
    field<vec> yF = List2Field_vec(y);
    List Xbetas = as<List>(Data["Xbetas"]);
    field<vec> XbetasF = List2Field_vec(Xbetas);
    List Z = as<List>(Data["Z"]);
    field<mat> ZF = List2Field_mat(Z);
    List RE_inds = as<List>(Data["RE_inds"]);
    field<uvec> RE_indsF = List2Field_uvec(RE_inds);
    List RE_inds2 = as<List>(Data["RE_inds2"]);
    field<uvec> RE_inds2F = List2Field_uvec(RE_inds2);
    List sigmas = as<List>(Data["sigmas"]);
    mat invD = as<mat>(Data["invD"]);
    List idL = as<List>(Data["idL"]);
    field<uvec> idLF = List2Field_uvec(idL);
    List idL2 = as<List>(Data["idL2"]);
    field<uvec> idL2F = List2Field_uvec(idL2, false);
    CharacterVector fams = as< CharacterVector>(Data["fams"]);
    CharacterVector links = as< CharacterVector>(Data["links"]);
    vec event = as<vec>(Data["event"]);
    uvec idGK = as<uvec>(Data["idGK_fast"]);
    mat W1 = as<mat>(Data["W1"]);
    mat W1s = as<mat>(Data["W1s"]);
    mat W2 = as<mat>(Data["W2"]);
    mat W2s = as<mat>(Data["W2s"]);
    List U = as<List>(Data["U"]);
    field<mat> UF = List2Field_mat(U);
    List Us = as<List>(Data["Us"]);
    field<mat> UsF = List2Field_mat(Us);
    uvec row_inds_U = as<uvec>(Data["row_inds_U"]) - 1;
    uvec row_inds_Us = as<uvec>(Data["row_inds_Us"]) - 1;
    List col_inds = as<List>(Data["col_inds"]);
    field<uvec> col_indsF = List2Field_uvec(col_inds);
    List XXbetas = as<List>(Data["XXbetas"]);
    field<vec> XXbetasF = List2Field_vec(XXbetas);
    List XXsbetas = as<List>(Data["XXsbetas"]);
    field<vec> XXsbetasF = List2Field_vec(XXsbetas);
    List ZZ = as<List>(Data["ZZ"]);
    field<mat> ZZF = List2Field_mat(ZZ);
    List ZZs = as<List>(Data["ZZs"]);
    field<mat> ZZsF = List2Field_mat(ZZs);
    List idT = as<List>(Data["idT"]);
    field<uvec> idTF = List2Field_uvec(idT);
    List idTs = as<List>(Data["idTs"]);
    field<uvec> idTsF = List2Field_uvec(idTs);
    vec Pw = as<vec>(Data["Pw"]);
    if (interval_cens) {
        mat W1s_int = as<mat>(Data["W1s_int"]);
        mat W2s_int = as<mat>(Data["W2s_int"]);
        List Us_int = as<List>(Data["Us_int"]);
        List XXsbetas_int = as<List>(Data["XXsbetas_int"]);
        field<vec> XXsbetasF_int = List2Field_vec(XXsbetas_int);
        field<mat> UsF_int = List2Field_mat(Us_int);
        List ZZs_int = as<List>(Data["ZZs_int"]);
        field<mat> ZZsF_int = List2Field_mat(ZZs_int);
        vec Pw_int = as<vec>(Data["Pw_int"]);
    }
    // Initials
    mat b = as<mat>(initials["b"]);
    vec Bs_gammas = as<vec>(initials["Bs_gammas"]);
    vec gammas = as<vec>(initials["gammas"]);
    vec alphas = as<vec>(initials["alphas"]);
    vec phi_gammas = as<vec>(initials["phi_gammas"]);
    vec phi_alphas = as<vec>(initials["phi_alphas"]);
    double tau_Bs_gammas = as<double>(initials["tau_Bs_gammas"]);
    double tau_gammas = as<double>(initials["tau_gammas"]);
    double tau_alphas = as<double>(initials["tau_alphas"]);
    int n = b.n_rows;
    int ns = Pw.n_rows;
    // Priors
    vec mean_Bs_gammas = as<vec>(priors["mean_Bs_gammas"]);
    mat Tau_Bs_gammas = as<mat>(priors["Tau_Bs_gammas"]);
    vec mean_gammas = as<vec>(priors["mean_gammas"]);
    mat Tau_gammas = as<mat>(priors["Tau_gammas"]);
    vec mean_alphas = as<vec>(priors["mean_alphas"]);
    mat Tau_alphas = as<mat>(priors["Tau_alphas"]);
    double A_tau_Bs_gammas = as<double>(priors["A_tau_Bs_gammas"]);
    double B_tau_Bs_gammas = as<double>(priors["B_tau_Bs_gammas"]);
    double Apost_tau_Bs_gammas = A_tau_Bs_gammas + 0.5 * as<double>(priors["rank_Tau_Bs_gammas"]);
    //
    bool shrink_gammas = as<bool>(priors["shrink_gammas"]);
    double A_tau_gammas = as<double>(priors["A_tau_gammas"]);
    double B_tau_gammas = as<double>(priors["B_tau_gammas"]);
    double Apost_tau_gammas = A_tau_gammas + 0.5 * as<double>(priors["rank_Tau_gammas"]);
    double A_phi_gammas = as<double>(priors["A_phi_gammas"]);
    double B_phi_gammas = as<double>(priors["B_phi_gammas"]);
    double Apost_phi_gammas = A_phi_gammas + 0.5;
    //
    bool shrink_alphas = as<bool>(priors["shrink_alphas"]);
    double A_tau_alphas = as<double>(priors["A_tau_alphas"]);
    double B_tau_alphas = as<double>(priors["B_tau_alphas"]);
    double Apost_tau_alphas = A_tau_alphas + 0.5 * as<double>(priors["rank_Tau_alphas"]);
    double A_phi_alphas = as<double>(priors["A_phi_alphas"]);
    double B_phi_alphas = as<double>(priors["B_phi_alphas"]);
    double Apost_phi_alphas = A_phi_alphas + 0.5;
    // Control parameters
    int n_iter = as<int>(control["n_iter"]);
    int n_burnin = as<int>(control["n_burnin"]);
    int n_block = as<int>(control["n_block"]);
    int n_thin = as<int>(control["n_thin"]);
    double target_acc = as<double>(control["target_acc"]);
    double eps1 = as<double>(control["eps1"]);
    double eps2 = as<double>(control["eps2"]);
    double eps3 = as<double>(control["eps3"]);
    double c0 = as<double>(control["c0"]);
    double mc1 = - as<double>(control["c1"]);
    bool adaptCov = as<bool>(control["adaptCov"]);
    // derived from control
    int total_it = n_iter + n_burnin;
    int its = ceil(total_it / n_block);
    int n_out = floor(n_iter / n_thin);
    int start_cov_update = floor(n_burnin / n_block) - 1;
    uvec keep = seqC(n_burnin + 1, total_it, n_thin);
    // initial scales
    int n_b = b.n_cols;
    int n_Bs_gammas = Bs_gammas.n_rows;
    int n_gammas = gammas.n_rows;
    int n_alphas = alphas.n_rows;
    vec sigma_b = as<vec>(scales["b"]);
    double sigma_Bs_gammas = as<double>(scales["Bs_gammas"]);
    double sigma_gammas = as<double>(scales["gammas"]);
    double sigma_alphas = as<double>(scales["alphas"]);
    cube Sigma_b = as<cube>(Covs["b"]);
    mat Sigma_Bs_gammas = as<mat>(Covs["Bs_gammas"]);
    mat Sigma_gammas = as<mat>(Covs["gammas"]);
    mat Sigma_alphas = as<mat>(Covs["alphas"]);
    // store results
    cube out_b = cube(n, n_b, n_out, fill::zeros);
    mat out_Bs_gammas(n_Bs_gammas, n_out, fill::zeros);
    mat out_gammas(n_gammas, n_out, fill::zeros);
    mat out_alphas(n_alphas, n_out, fill::zeros);
    vec out_tau_Bs_gammas(n_out, fill::zeros);
    vec out_tau_gammas(n_out, fill::zeros);
    vec out_tau_alphas(n_out, fill::zeros);
    mat out_phi_gammas(n_gammas, n_out, fill::zeros);
    mat out_phi_alphas(n_alphas, n_out, fill::zeros);
    vec out_logWeights(n_out, fill::zeros);
    RNGScope scope;
    int keep_iterator = 0;
    int check_iterator = 0;
    vec log_us = log(randu<vec>(its * n_block));
    mat log_us_RE = log(randu<mat>(n, its * n_block));
    for (int it = 0; it < its; ++it) {
        mat accept_b(n, n_block, fill::zeros);
        vec accept_s(n_block, fill::zeros);
        cube block_b(n, n_b, n_block, fill::zeros);
        mat block_Bs_gammas(n_Bs_gammas, n_block, fill::zeros);
        mat block_gammas(n_gammas, n_block, fill::zeros);
        mat block_alphas(n_alphas, n_block, fill::zeros);
        cube rand_b = mvrnorm_cube(n_block, Sigma_b, sigma_b);
        mat rand_Bs_gammas = mvrnorm(n_block, sigma_Bs_gammas * Sigma_Bs_gammas);
        mat rand_gammas = mvrnorm(n_block, sigma_gammas * Sigma_gammas);
        mat rand_alphas = mvrnorm(n_block, sigma_alphas * Sigma_alphas);
        mat new_b(n, n_b, fill::zeros);
        vec new_Bs_gammas(n_Bs_gammas, fill::zeros);
        vec new_gammas(n_gammas, fill::zeros);
        vec new_alphas(n_alphas, fill::zeros);
        for (int i = 0; i < n_block; ++i) {
            int iter = it * n_block + i;
            // sample random effects
            vec current_logpost_RE = log_postREF(b, Bs_gammas, gammas, alphas, yF, XbetasF, ZF,
                                                 RE_indsF, RE_inds2F, idLF, idL2F, fams,
                                                 links, sigmas, invD, n, ns, n_alphas,
                                                 event, W1, W1s, W2, W2s, XXbetasF,
                                                 XXsbetasF, ZZF, ZZsF, UF, UsF,
                                                 Pw, idGK, idTF, idTsF, col_indsF,
                                                 row_inds_U, row_inds_Us);

            new_b = b + extract_b(rand_b, i);
            vec new_logpost_RE = log_postREF(new_b, Bs_gammas, gammas, alphas, yF, XbetasF, ZF,
                                             RE_indsF, RE_inds2F, idLF, idL2F, fams,
                                             links, sigmas, invD, n, ns, n_alphas,
                                             event, W1, W1s, W2, W2s, XXbetasF,
                                             XXsbetasF, ZZF, ZZsF, UF, UsF,
                                             Pw, idGK, idTF, idTsF, col_indsF,
                                             row_inds_U, row_inds_Us);
            vec lRatio_RE = new_logpost_RE - current_logpost_RE;
            for (int m = 0; m < n; ++m) {
                if (log_us_RE.at(m, iter) < lRatio_RE[m]) {
                    accept_b.at(m, i) = 1.0;
                    b.row(m) = new_b.row(m);
                }
            }
            // sample survival
            mat Wlong = lin_pred_matF(XXbetasF, ZZF, b, UF, RE_inds2F, idTF,
                                      col_indsF, row_inds_U, n, n_alphas);
            mat Wlongs = lin_pred_matF(XXsbetasF, ZZsF, b, UsF, RE_inds2F, idTsF,
                                       col_indsF, row_inds_Us, ns, n_alphas);
            double current_logpost = logPosterior(event, W1, W1s, Bs_gammas,
                                                  W2, W2s, gammas, Wlong, Wlongs, alphas, Pw,
                                                  mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                                  mean_gammas, Tau_gammas, tau_gammas,
                                                  mean_alphas, Tau_alphas, tau_alphas);
            new_Bs_gammas = Bs_gammas + rand_Bs_gammas.col(i);
            new_gammas = gammas + rand_gammas.col(i);
            new_alphas = alphas + rand_alphas.col(i);
            double new_logpost = logPosterior(event, W1, W1s, new_Bs_gammas,
                                              W2, W2s, new_gammas, Wlong, Wlongs, new_alphas, Pw,
                                              mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                              mean_gammas, Tau_gammas, tau_gammas,
                                              mean_alphas, Tau_alphas, tau_alphas);
            double lRatio = new_logpost - current_logpost;
            if (log_us.at(iter) < lRatio) {
                accept_s.at(i) = 1.0;
                Bs_gammas = new_Bs_gammas;
                gammas = new_gammas;
                alphas = new_alphas;
            }
            block_b.slice(i) = b;
            block_Bs_gammas.col(i) = Bs_gammas;
            block_gammas.col(i) = gammas;
            block_alphas.col(i) = alphas;
            // sample scales
            double Bpost_tau_Bs_gammas = B_tau_Bs_gammas + 0.5 * as_scalar(Bs_gammas.t() * Tau_Bs_gammas * Bs_gammas);
            tau_Bs_gammas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_Bs_gammas, 1.0 / Bpost_tau_Bs_gammas)));
            if (shrink_gammas) {
                for (int k1 = 0; k1 < n_gammas; ++k1) {
                    double Bpost_phi_gammas = B_phi_gammas + 0.5 * tau_gammas * pow(gammas.at(k1), 2);
                    phi_gammas.at(k1) = std::min(eps3,
                                     std::max(eps2, ::Rf_rgamma(Apost_phi_gammas, 1.0 / Bpost_phi_gammas)));
                }
                Tau_gammas.diag() = phi_gammas;
                double Bpost_tau_gammas = B_tau_gammas + 0.5 * as_scalar(gammas.t() * Tau_gammas * gammas);
                tau_gammas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_gammas, 1.0 / Bpost_tau_gammas)));
            }
            if (shrink_alphas) {
                for (int k2 = 0; k2 < n_alphas; ++k2) {
                    double Bpost_phi_alphas = B_phi_alphas + 0.5 * tau_alphas * pow(alphas.at(k2), 2);
                    phi_alphas.at(k2) = std::min(eps3,
                                  std::max(eps2, ::Rf_rgamma(Apost_phi_alphas, 1.0 / Bpost_phi_alphas)));
                }
                Tau_alphas.diag() = phi_alphas;
                double Bpost_tau_alphas = B_tau_alphas + 0.5 * as_scalar(alphas.t() * Tau_alphas * alphas);
                tau_alphas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_alphas, 1.0 / Bpost_tau_alphas)));
            }
            // store results
            if ((unsigned)iter == keep[check_iterator]) {
                out_b.slice(keep_iterator) = b;
                out_Bs_gammas.col(keep_iterator) = Bs_gammas;
                out_gammas.col(keep_iterator) = gammas;
                out_alphas.col(keep_iterator) = alphas;
                out_phi_gammas.col(keep_iterator) = phi_gammas;
                out_phi_alphas.col(keep_iterator) = phi_alphas;
                out_tau_Bs_gammas.at(keep_iterator) = tau_Bs_gammas;
                out_tau_gammas.at(keep_iterator) = tau_gammas;
                out_tau_alphas.at(keep_iterator) = tau_alphas;
                out_logWeights.at(keep_iterator) = log_weightsREF(b, yF, XbetasF, ZF, RE_indsF,
                                  idLF, idL2F, fams, links, sigmas, invD, n);
                keep_iterator += 1;
                check_iterator += 1;
            }
        }
        vec rate_b = mean(accept_b, 1);
        double rate_s = mean(accept_s);
        double g1 = pow(it + 1, mc1);
        double g2 = c0 * g1;
        for (int m = 0; m < n; ++m) {
            sigma_b[m] = bounds_sigma(exp(log(sigma_b[m]) + g2 * (rate_b[m] - target_acc)), eps1, eps3);
        }
        sigma_Bs_gammas = bounds_sigma(exp(log(sigma_Bs_gammas) + g2 * (rate_s - target_acc)), eps1, eps3);
        sigma_gammas = bounds_sigma(exp(log(sigma_gammas) + g2 * (rate_s - target_acc)), eps1, eps3);
        sigma_alphas = bounds_sigma(exp(log(sigma_alphas) + g2 * (rate_s - target_acc)), eps1, eps3);
        if (adaptCov && it > start_cov_update) {
            Sigma_Bs_gammas = bounds_Cov(Sigma_Bs_gammas + g1 * (cov(block_Bs_gammas.t()) - Sigma_Bs_gammas), eps2, eps3);
            Sigma_gammas = bounds_Cov(Sigma_gammas + g1 * (cov(block_gammas.t()) - Sigma_gammas), eps2, eps3);
            Sigma_alphas = bounds_Cov(Sigma_alphas + g1 * (cov(block_alphas.t()) - Sigma_alphas), eps2, eps3);
        }
    }
    // export
    return List::create(Named("mcmc") = List::create(
        Named("b") = out_b,
        Named("Bs_gammas") = out_Bs_gammas,
        Named("gammas") = out_gammas,
        Named("alphas") = out_alphas,
        Named("phi_gammas") = out_phi_gammas,
        Named("phi_alphas") = out_phi_alphas,
        Named("tau_Bs_gammas") = out_tau_Bs_gammas,
        Named("tau_gammas") = out_tau_gammas,
        Named("tau_alphas") = out_tau_alphas
    ),
    Named("logWeights") = out_logWeights,
    Named("scales") = List::create(
        Named("sigma") = List::create(
            Named("b") = sigma_b,
            Named("Bs_gammas") = sigma_Bs_gammas,
            Named("gammas") = sigma_gammas,
            Named("alphas") = sigma_alphas
        ),
        Named("Sigma") = List::create(
            Named("Bs_gammas") = Sigma_Bs_gammas,
            Named("gammas") = Sigma_gammas,
            Named("alphas") = Sigma_alphas
        )
    ));
}

// [[Rcpp::export]]
List lap_rwm_C_nogammas (List initials, List Data, List priors, List scales, List Covs,
                         List control) {
    // Data
    List y = as<List>(Data["y"]);
    field<vec> yF = List2Field_vec(y);
    List Xbetas = as<List>(Data["Xbetas"]);
    field<vec> XbetasF = List2Field_vec(Xbetas);
    List Z = as<List>(Data["Z"]);
    field<mat> ZF = List2Field_mat(Z);
    List RE_inds = as<List>(Data["RE_inds"]);
    field<uvec> RE_indsF = List2Field_uvec(RE_inds);
    List RE_inds2 = as<List>(Data["RE_inds2"]);
    field<uvec> RE_inds2F = List2Field_uvec(RE_inds2);
    List sigmas = as<List>(Data["sigmas"]);
    mat invD = as<mat>(Data["invD"]);
    List idL = as<List>(Data["idL"]);
    field<uvec> idLF = List2Field_uvec(idL);
    List idL2 = as<List>(Data["idL2"]);
    field<uvec> idL2F = List2Field_uvec(idL2, false);
    CharacterVector fams = as< CharacterVector>(Data["fams"]);
    CharacterVector links = as< CharacterVector>(Data["links"]);
    vec event = as<vec>(Data["event"]);
    uvec idGK = as<uvec>(Data["idGK_fast"]);
    mat W1 = as<mat>(Data["W1"]);
    mat W1s = as<mat>(Data["W1s"]);
    List U = as<List>(Data["U"]);
    field<mat> UF = List2Field_mat(U);
    List Us = as<List>(Data["Us"]);
    field<mat> UsF = List2Field_mat(Us);
    uvec row_inds_U = as<uvec>(Data["row_inds_U"]) - 1;
    uvec row_inds_Us = as<uvec>(Data["row_inds_Us"]) - 1;
    List col_inds = as<List>(Data["col_inds"]);
    field<uvec> col_indsF = List2Field_uvec(col_inds);
    List XXbetas = as<List>(Data["XXbetas"]);
    field<vec> XXbetasF = List2Field_vec(XXbetas);
    List XXsbetas = as<List>(Data["XXsbetas"]);
    field<vec> XXsbetasF = List2Field_vec(XXsbetas);
    List ZZ = as<List>(Data["ZZ"]);
    field<mat> ZZF = List2Field_mat(ZZ);
    List ZZs = as<List>(Data["ZZs"]);
    field<mat> ZZsF = List2Field_mat(ZZs);
    List idT = as<List>(Data["idT"]);
    field<uvec> idTF = List2Field_uvec(idT);
    List idTs = as<List>(Data["idTs"]);
    field<uvec> idTsF = List2Field_uvec(idTs);
    vec Pw = as<vec>(Data["Pw"]);
    // Initials
    mat b = as<mat>(initials["b"]);
    vec Bs_gammas = as<vec>(initials["Bs_gammas"]);
    vec alphas = as<vec>(initials["alphas"]);
    vec phi_alphas = as<vec>(initials["phi_alphas"]);
    double tau_Bs_gammas = as<double>(initials["tau_Bs_gammas"]);
    double tau_alphas = as<double>(initials["tau_alphas"]);
    int n = b.n_rows;
    int ns = Pw.n_rows;
    // Priors
    vec mean_Bs_gammas = as<vec>(priors["mean_Bs_gammas"]);
    mat Tau_Bs_gammas = as<mat>(priors["Tau_Bs_gammas"]);
    vec mean_alphas = as<vec>(priors["mean_alphas"]);
    mat Tau_alphas = as<mat>(priors["Tau_alphas"]);
    double A_tau_Bs_gammas = as<double>(priors["A_tau_Bs_gammas"]);
    double B_tau_Bs_gammas = as<double>(priors["B_tau_Bs_gammas"]);
    double Apost_tau_Bs_gammas = A_tau_Bs_gammas + 0.5 * as<double>(priors["rank_Tau_Bs_gammas"]);
    //
    bool shrink_alphas = as<bool>(priors["shrink_alphas"]);
    double A_tau_alphas = as<double>(priors["A_tau_alphas"]);
    double B_tau_alphas = as<double>(priors["B_tau_alphas"]);
    double Apost_tau_alphas = A_tau_alphas + 0.5 * as<double>(priors["rank_Tau_alphas"]);
    double A_phi_alphas = as<double>(priors["A_phi_alphas"]);
    double B_phi_alphas = as<double>(priors["B_phi_alphas"]);
    double Apost_phi_alphas = A_phi_alphas + 0.5;
    // Control parameters
    int n_iter = as<int>(control["n_iter"]);
    int n_burnin = as<int>(control["n_burnin"]);
    int n_block = as<int>(control["n_block"]);
    int n_thin = as<int>(control["n_thin"]);
    double target_acc = as<double>(control["target_acc"]);
    double eps1 = as<double>(control["eps1"]);
    double eps2 = as<double>(control["eps2"]);
    double eps3 = as<double>(control["eps3"]);
    double c0 = as<double>(control["c0"]);
    double mc1 = - as<double>(control["c1"]);
    bool adaptCov = as<bool>(control["adaptCov"]);
    // derived from control
    int total_it = n_iter + n_burnin;
    int its = ceil(total_it / n_block);
    int n_out = floor(n_iter / n_thin);
    int start_cov_update = floor(n_burnin / n_block) - 1;
    uvec keep = seqC(n_burnin + 1, total_it, n_thin);
    // initial scales
    int n_b = b.n_cols;
    int n_Bs_gammas = Bs_gammas.n_rows;
    int n_alphas = alphas.n_rows;
    vec sigma_b = as<vec>(scales["b"]);
    double sigma_Bs_gammas = as<double>(scales["Bs_gammas"]);
    double sigma_alphas = as<double>(scales["alphas"]);
    cube Sigma_b = as<cube>(Covs["b"]);
    mat Sigma_Bs_gammas = as<mat>(Covs["Bs_gammas"]);
    mat Sigma_alphas = as<mat>(Covs["alphas"]);
    // store results
    cube out_b = cube(n, n_b, n_out, fill::zeros);
    mat out_Bs_gammas(n_Bs_gammas, n_out, fill::zeros);
    mat out_alphas(n_alphas, n_out, fill::zeros);
    vec out_tau_Bs_gammas(n_out, fill::zeros);
    vec out_tau_alphas(n_out, fill::zeros);
    mat out_phi_alphas(n_alphas, n_out, fill::zeros);
    vec out_logWeights(n_out, fill::zeros);
    RNGScope scope;
    int keep_iterator = 0;
    int check_iterator = 0;
    vec log_us = log(randu<vec>(its * n_block));
    mat log_us_RE = log(randu<mat>(n, its * n_block));
    for (int it = 0; it < its; ++it) {
        mat accept_b(n, n_block, fill::zeros);
        vec accept_s(n_block, fill::zeros);
        cube block_b(n, n_b, n_block, fill::zeros);
        mat block_Bs_gammas(n_Bs_gammas, n_block, fill::zeros);
        mat block_alphas(n_alphas, n_block, fill::zeros);
        cube rand_b = mvrnorm_cube(n_block, Sigma_b, sigma_b);
        mat rand_Bs_gammas = mvrnorm(n_block, sigma_Bs_gammas * Sigma_Bs_gammas);
        mat rand_alphas = mvrnorm(n_block, sigma_alphas * Sigma_alphas);
        mat new_b(n, n_b, fill::zeros);
        vec new_Bs_gammas(n_Bs_gammas, fill::zeros);
        vec new_alphas(n_alphas, fill::zeros);
        for (int i = 0; i < n_block; ++i) {
            int iter = it * n_block + i;
            // sample random effects
            vec current_logpost_RE = log_postRE_nogammasF(b, Bs_gammas, alphas, yF, XbetasF, ZF,
                                                          RE_indsF, RE_inds2F, idLF, idL2F, fams,
                                                          links, sigmas, invD, n, ns, n_alphas, event, W1, W1s,
                                                          XXbetasF, XXsbetasF, ZZF, ZZsF, UF, UsF,
                                                          Pw, idGK, idTF, idTsF, col_indsF,
                                                          row_inds_U, row_inds_Us);

            new_b = b + extract_b(rand_b, i);
            vec new_logpost_RE = log_postRE_nogammasF(new_b, Bs_gammas, alphas, yF, XbetasF, ZF,
                                                      RE_indsF, RE_inds2F, idLF, idL2F, fams,
                                                      links, sigmas, invD, n, ns, n_alphas, event, W1, W1s,
                                                      XXbetasF, XXsbetasF, ZZF, ZZsF, UF, UsF,
                                                      Pw, idGK, idTF, idTsF, col_indsF,
                                                      row_inds_U, row_inds_Us);
            vec lRatio_RE = new_logpost_RE - current_logpost_RE;
            for (int m = 0; m < n; ++m) {
                if (log_us_RE.at(m, iter) < lRatio_RE[m]) {
                    accept_b.at(m, i) = 1.0;
                    b.row(m) = new_b.row(m);
                }
            }
            // sample survival
            mat Wlong = lin_pred_matF(XXbetasF, ZZF, b, UF, RE_inds2F, idTF,
                                      col_indsF, row_inds_U, n, n_alphas);
            mat Wlongs = lin_pred_matF(XXsbetasF, ZZsF, b, UsF, RE_inds2F, idTsF,
                                       col_indsF, row_inds_Us, ns, n_alphas);
            double current_logpost = logPosterior_nogammas(event, W1, W1s, Bs_gammas,
                                                           Wlong, Wlongs, alphas, Pw,
                                                           mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                                           mean_alphas, Tau_alphas, tau_alphas);
            new_Bs_gammas = Bs_gammas + rand_Bs_gammas.col(i);
            new_alphas = alphas + rand_alphas.col(i);
            double new_logpost = logPosterior_nogammas(event, W1, W1s, new_Bs_gammas,
                                                       Wlong, Wlongs, new_alphas, Pw,
                                                       mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                                       mean_alphas, Tau_alphas, tau_alphas);
            double lRatio = new_logpost - current_logpost;
            if (log_us.at(iter) < lRatio) {
                accept_s.at(i) = 1.0;
                Bs_gammas = new_Bs_gammas;
                alphas = new_alphas;
            }
            block_b.slice(i) = b;
            block_Bs_gammas.col(i) = Bs_gammas;
            block_alphas.col(i) = alphas;
            // sample scales
            double Bpost_tau_Bs_gammas = B_tau_Bs_gammas + 0.5 * as_scalar(Bs_gammas.t() * Tau_Bs_gammas * Bs_gammas);
            tau_Bs_gammas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_Bs_gammas, 1.0 / Bpost_tau_Bs_gammas)));
            if (shrink_alphas) {
                for (int k2 = 0; k2 < n_alphas; ++k2) {
                    double Bpost_phi_alphas = B_phi_alphas + 0.5 * tau_alphas * pow(alphas.at(k2), 2);
                    phi_alphas.at(k2) = std::min(eps3,
                                  std::max(eps2, ::Rf_rgamma(Apost_phi_alphas, 1.0 / Bpost_phi_alphas)));
                }
                Tau_alphas.diag() = phi_alphas;
                double Bpost_tau_alphas = B_tau_alphas + 0.5 * as_scalar(alphas.t() * Tau_alphas * alphas);
                tau_alphas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_alphas, 1.0 / Bpost_tau_alphas)));
            }
            // store results
            if ((unsigned)iter == keep[check_iterator]) {
                out_b.slice(keep_iterator) = b;
                out_Bs_gammas.col(keep_iterator) = Bs_gammas;
                out_alphas.col(keep_iterator) = alphas;
                out_phi_alphas.col(keep_iterator) = phi_alphas;
                out_tau_Bs_gammas.at(keep_iterator) = tau_Bs_gammas;
                out_tau_alphas.at(keep_iterator) = tau_alphas;
                out_logWeights.at(keep_iterator) = log_weightsREF(b, yF, XbetasF, ZF, RE_indsF,
                                  idLF, idL2F, fams, links, sigmas, invD, n);
                keep_iterator += 1;
                check_iterator += 1;
            }
        }
        vec rate_b = mean(accept_b, 1);
        double rate_s = mean(accept_s);
        double g1 = pow(it + 1, mc1);
        double g2 = c0 * g1;
        for (int m = 0; m < n; ++m) {
            sigma_b[m] = bounds_sigma(exp(log(sigma_b[m]) + g2 * (rate_b[m] - target_acc)), eps1, eps3);
        }
        sigma_Bs_gammas = bounds_sigma(exp(log(sigma_Bs_gammas) + g2 * (rate_s - target_acc)), eps1, eps3);
        sigma_alphas = bounds_sigma(exp(log(sigma_alphas) + g2 * (rate_s - target_acc)), eps1, eps3);
        if (adaptCov && it > start_cov_update) {
            Sigma_Bs_gammas = bounds_Cov(Sigma_Bs_gammas + g1 * (cov(block_Bs_gammas.t()) - Sigma_Bs_gammas), eps2, eps3);
            Sigma_alphas = bounds_Cov(Sigma_alphas + g1 * (cov(block_alphas.t()) - Sigma_alphas), eps2, eps3);
        }
    }
    // export
    return List::create(Named("mcmc") = List::create(
        Named("b") = out_b,
        Named("Bs_gammas") = out_Bs_gammas,
        Named("alphas") = out_alphas,
        Named("phi_alphas") = out_phi_alphas,
        Named("tau_Bs_gammas") = out_tau_Bs_gammas,
        Named("tau_alphas") = out_tau_alphas
    ),
    Named("logWeights") = out_logWeights,
    Named("scales") = List::create(
        Named("sigma") = List::create(
            Named("b") = sigma_b,
            Named("Bs_gammas") = sigma_Bs_gammas,
            Named("alphas") = sigma_alphas
        ),
        Named("Sigma") = List::create(
            Named("Bs_gammas") = Sigma_Bs_gammas,
            Named("alphas") = Sigma_alphas
        )
    ));
}

double logPosterior_woRE (const vec& event,
                          const mat& W1, const mat& W1s, const vec& Bs_gammas,
                          const mat& W2, const mat& W2s, const vec& gammas,
                          const mat& Wlong, const mat& Wlongs, const vec& alphas,
                          const vec& Pw,
                          const vec& mean_Bs_gammas, const mat& Tau_Bs_gammas, double tau_Bs_gammas,
                          const vec& mean_gammas, const mat& Tau_gammas, double tau_gammas,
                          const vec& mean_alphas, const mat& Tau_alphas, double tau_alphas) {
    // log-likelihood
    vec log_h = W1 * Bs_gammas + W2 * gammas + Wlong * alphas;
    vec H = Pw % exp(W1s * Bs_gammas + W2s * gammas + Wlongs * alphas);
    //
    double out = sum(event % log_h) - sum(H) +
        logPrior(Bs_gammas, mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas) +
        logPrior(gammas, mean_gammas, Tau_gammas, tau_gammas) +
        logPrior(alphas, mean_alphas, Tau_alphas, tau_alphas);
    return(out);
}

double logPosterior_woRE_nogammas (const vec& event,
                                   const mat& W1, const mat& W1s, const vec& Bs_gammas,
                                   const mat& Wlong, const mat& Wlongs, const vec& alphas,
                                   const vec& Pw,
                                   const vec& mean_Bs_gammas, const mat& Tau_Bs_gammas, double tau_Bs_gammas,
                                   const vec& mean_alphas, const mat& Tau_alphas, double tau_alphas) {
    // log-likelihood
    vec log_h = W1 * Bs_gammas + Wlong * alphas;
    vec H = Pw % exp(W1s * Bs_gammas + Wlongs * alphas);
    //
    double out = sum(event % log_h) - sum(H) +
        logPrior(Bs_gammas, mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas) +
        logPrior(alphas, mean_alphas, Tau_alphas, tau_alphas);
    return(out);
}

// [[Rcpp::export]]
List lap_rwm_C_woRE (List initials, List Data, List priors, List scales, List Covs, 
                     List control) {
    // Data
    vec event = as<vec>(Data["event"]);
    uvec idGK_fast = as<uvec>(Data["idGK_fast"]);
    mat W1 = as<mat>(Data["W1"]);
    mat W1s = as<mat>(Data["W1s"]);
    mat W2 = as<mat>(Data["W2"]);
    mat W2s = as<mat>(Data["W2s"]);
    mat Wlong = as<mat>(Data["Wlong"]);
    mat Wlongs = as<mat>(Data["Wlongs"]);
    vec Pw = as<vec>(Data["Pw"]);
    // Initials
    vec Bs_gammas = as<vec>(initials["Bs_gammas"]);
    vec gammas = as<vec>(initials["gammas"]);
    vec alphas = as<vec>(initials["alphas"]);
    vec phi_gammas = as<vec>(initials["phi_gammas"]);
    vec phi_alphas = as<vec>(initials["phi_alphas"]);
    double tau_Bs_gammas = as<double>(initials["tau_Bs_gammas"]);
    double tau_gammas = as<double>(initials["tau_gammas"]);
    double tau_alphas = as<double>(initials["tau_alphas"]);
    // Priors
    vec mean_Bs_gammas = as<vec>(priors["mean_Bs_gammas"]);
    mat Tau_Bs_gammas = as<mat>(priors["Tau_Bs_gammas"]);
    vec mean_gammas = as<vec>(priors["mean_gammas"]);
    mat Tau_gammas = as<mat>(priors["Tau_gammas"]);
    vec mean_alphas = as<vec>(priors["mean_alphas"]);
    mat Tau_alphas = as<mat>(priors["Tau_alphas"]);
    double A_tau_Bs_gammas = as<double>(priors["A_tau_Bs_gammas"]);
    double B_tau_Bs_gammas = as<double>(priors["B_tau_Bs_gammas"]);
    double Apost_tau_Bs_gammas = A_tau_Bs_gammas + 0.5 * as<double>(priors["rank_Tau_Bs_gammas"]);
    //
    bool shrink_gammas = as<bool>(priors["shrink_gammas"]);
    double A_tau_gammas = as<double>(priors["A_tau_gammas"]);
    double B_tau_gammas = as<double>(priors["B_tau_gammas"]);
    double Apost_tau_gammas = A_tau_gammas + 0.5 * as<double>(priors["rank_Tau_gammas"]);
    double A_phi_gammas = as<double>(priors["A_phi_gammas"]);
    double B_phi_gammas = as<double>(priors["B_phi_gammas"]);
    double Apost_phi_gammas = A_phi_gammas + 0.5;
    //
    bool shrink_alphas = as<bool>(priors["shrink_alphas"]);
    double A_tau_alphas = as<double>(priors["A_tau_alphas"]);
    double B_tau_alphas = as<double>(priors["B_tau_alphas"]);
    double Apost_tau_alphas = A_tau_alphas + 0.5 * as<double>(priors["rank_Tau_alphas"]);
    double A_phi_alphas = as<double>(priors["A_phi_alphas"]);
    double B_phi_alphas = as<double>(priors["B_phi_alphas"]);
    double Apost_phi_alphas = A_phi_alphas + 0.5;
    // Control parameters
    int n_iter = as<int>(control["n_iter"]);
    int n_burnin = as<int>(control["n_burnin"]);
    int n_block = as<int>(control["n_block"]);
    int n_thin = as<int>(control["n_thin"]);
    double target_acc = as<double>(control["target_acc"]);
    double eps1 = as<double>(control["eps1"]);
    double eps2 = as<double>(control["eps2"]);
    double eps3 = as<double>(control["eps3"]);
    double c0 = as<double>(control["c0"]);
    double mc1 = - as<double>(control["c1"]);
    bool adaptCov = as<bool>(control["adaptCov"]);
    // derived from control
    int total_it = n_iter + n_burnin;
    int its = ceil(total_it / n_block);
    int n_out = floor(n_iter / n_thin);
    int start_cov_update = floor(n_burnin / n_block) - 1;
    uvec keep = seqC(n_burnin + 1, total_it, n_thin);
    // initial
    int n_Bs_gammas = Bs_gammas.n_rows;
    int n_gammas = gammas.n_rows;
    int n_alphas = alphas.n_rows;
    double sigma_Bs_gammas = as<double>(scales["Bs_gammas"]);
    double sigma_gammas = as<double>(scales["gammas"]);
    double sigma_alphas = as<double>(scales["alphas"]);
    mat Sigma_Bs_gammas = as<mat>(Covs["Bs_gammas"]);
    mat Sigma_gammas = as<mat>(Covs["gammas"]);
    mat Sigma_alphas = as<mat>(Covs["alphas"]);
    // store results
    mat out_Bs_gammas(n_Bs_gammas, n_out, fill::zeros);
    mat out_gammas(n_gammas, n_out, fill::zeros);
    mat out_alphas(n_alphas, n_out, fill::zeros);
    vec out_tau_Bs_gammas(n_out, fill::zeros);
    vec out_tau_gammas(n_out, fill::zeros);
    vec out_tau_alphas(n_out, fill::zeros);
    mat out_phi_gammas(n_gammas, n_out, fill::zeros);
    mat out_phi_alphas(n_alphas, n_out, fill::zeros);
    vec out_logPost(n_out, fill::zeros);
    RNGScope scope;
    int keep_iterator = 0;
    int check_iterator = 0;
    vec log_us = log(randu<vec>(its * n_block));
    double current_logpost = logPosterior_woRE(event, W1, W1s, Bs_gammas,
                                               W2, W2s, gammas, Wlong, Wlongs, alphas, Pw,
                                               mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                               mean_gammas, Tau_gammas, tau_gammas,
                                               mean_alphas, Tau_alphas, tau_alphas);
    for (int it = 0; it < its; ++it) {
        vec accept = vec(n_block, fill::zeros);
        mat block_Bs_gammas = mat(n_Bs_gammas, n_block, fill::zeros);
        mat block_gammas = mat(n_gammas, n_block, fill::zeros);
        mat block_alphas = mat(n_alphas, n_block, fill::zeros);
        mat rand_Bs_gammas = mvrnorm(n_block, sigma_Bs_gammas * Sigma_Bs_gammas);
        mat rand_gammas = mvrnorm(n_block, sigma_gammas * Sigma_gammas);
        mat rand_alphas = mvrnorm(n_block, sigma_alphas * Sigma_alphas);
        vec new_Bs_gammas = vec(n_Bs_gammas, fill::zeros);
        vec new_gammas = vec(n_gammas, fill::zeros);
        vec new_alphas = vec(n_alphas, fill::zeros);
        for (int i = 0; i < n_block; ++i) {
            new_Bs_gammas = Bs_gammas + rand_Bs_gammas.col(i);
            new_gammas = gammas + rand_gammas.col(i);
            new_alphas = alphas + rand_alphas.col(i);
            double new_logpost = logPosterior_woRE(event, W1, W1s, new_Bs_gammas,
                                                   W2, W2s, new_gammas, Wlong, Wlongs, 
                                                   new_alphas, Pw,
                                                   mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                                   mean_gammas, Tau_gammas, tau_gammas,
                                                   mean_alphas, Tau_alphas, tau_alphas);
            double lRatio = new_logpost - current_logpost;
            int iter = it * n_block + i;
            if (log_us.at(iter) < lRatio) {
                accept.at(i) = 1;
                Bs_gammas = new_Bs_gammas;
                gammas = new_gammas;
                alphas = new_alphas;
                current_logpost = new_logpost;
            }
            block_Bs_gammas.col(i) = Bs_gammas;
            block_gammas.col(i) = gammas;
            block_alphas.col(i) = alphas;
            // sample scales
            double Bpost_tau_Bs_gammas = B_tau_Bs_gammas + 0.5 * as_scalar(Bs_gammas.t() * Tau_Bs_gammas * Bs_gammas);
            tau_Bs_gammas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_Bs_gammas, 1.0 / Bpost_tau_Bs_gammas)));
            if (shrink_gammas) {
                for (int k1 = 0; k1 < n_gammas; ++k1) {
                    double Bpost_phi_gammas = B_phi_gammas + 0.5 * tau_gammas * pow(gammas.at(k1), 2);
                    phi_gammas.at(k1) = std::min(eps3,
                                  std::max(eps2, ::Rf_rgamma(Apost_phi_gammas, 1.0 / Bpost_phi_gammas)));
                }
                Tau_gammas.diag() = phi_gammas;
                double Bpost_tau_gammas = B_tau_gammas + 0.5 * as_scalar(gammas.t() * Tau_gammas * gammas);
                tau_gammas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_gammas, 1.0 / Bpost_tau_gammas)));
            }
            if (shrink_alphas) {
                for (int k2 = 0; k2 < n_alphas; ++k2) {
                    double Bpost_phi_alphas = B_phi_alphas + 0.5 * tau_alphas * pow(alphas.at(k2), 2);
                    phi_alphas.at(k2) = std::min(eps3,
                                  std::max(eps2, ::Rf_rgamma(Apost_phi_alphas, 1.0 / Bpost_phi_alphas)));
                }
                Tau_alphas.diag() = phi_alphas;
                double Bpost_tau_alphas = B_tau_alphas + 0.5 * as_scalar(alphas.t() * Tau_alphas * alphas);
                tau_alphas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_alphas, 1.0 / Bpost_tau_alphas)));
            }
            // store results
            if ((unsigned)iter == keep[check_iterator]) {
                out_Bs_gammas.col(keep_iterator) = Bs_gammas;
                out_gammas.col(keep_iterator) = gammas;
                out_alphas.col(keep_iterator) = alphas;
                out_phi_gammas.col(keep_iterator) = phi_gammas;
                out_phi_alphas.col(keep_iterator) = phi_alphas;
                out_tau_Bs_gammas.at(keep_iterator) = tau_Bs_gammas;
                out_tau_gammas.at(keep_iterator) = tau_gammas;
                out_tau_alphas.at(keep_iterator) = tau_alphas;
                out_logPost.at(keep_iterator) = current_logpost;
                keep_iterator += 1;
                check_iterator += 1;
            }
        }
        double rate = mean(accept);
        double g1 = pow(it + 1, mc1);
        double g2 = c0 * g1;
        sigma_Bs_gammas = bounds_sigma(exp(log(sigma_Bs_gammas) + g2 * (rate - target_acc)), eps1, eps3);
        sigma_gammas = bounds_sigma(exp(log(sigma_gammas) + g2 * (rate - target_acc)), eps1, eps3);
        sigma_alphas = bounds_sigma(exp(log(sigma_alphas) + g2 * (rate - target_acc)), eps1, eps3);
        if (adaptCov && it > start_cov_update) {
            Sigma_Bs_gammas = bounds_Cov(Sigma_Bs_gammas + g1 * (cov(block_Bs_gammas.t()) - Sigma_Bs_gammas), eps2, eps3);
            Sigma_gammas = bounds_Cov(Sigma_gammas + g1 * (cov(block_gammas.t()) - Sigma_gammas), eps2, eps3);
            Sigma_alphas = bounds_Cov(Sigma_alphas + g1 * (cov(block_alphas.t()) - Sigma_alphas), eps2, eps3);
        }
    }
    // export
    return List::create(
        Named("mcmc") = List::create(
            Named("Bs_gammas") = out_Bs_gammas,
            Named("gammas") = out_gammas,
            Named("alphas") = out_alphas,
            Named("phi_gammas") = out_phi_gammas,
            Named("phi_alphas") = out_phi_alphas,
            Named("tau_Bs_gammas") = out_tau_Bs_gammas,
            Named("tau_gammas") = out_tau_gammas,
            Named("tau_alphas") = out_tau_alphas
        ),
        
        Named("logPost") = out_logPost,
        
        Named("scales") = List::create(
            Named("sigma") = List::create(
                Named("Bs_gammas") = sigma_Bs_gammas,
                Named("gammas") = sigma_gammas,
                Named("alphas") = sigma_alphas
            ),
            Named("Sigma") = List::create(
                Named("Bs_gammas") = Sigma_Bs_gammas,
                Named("gammas") = Sigma_gammas,
                Named("alphas") = Sigma_alphas
            )
        )
    );
}

// [[Rcpp::export]]
List lap_rwm_C_woRE_nogammas (List initials, List Data, List priors, List scales, 
                              List Covs, List control) {
    // Data
    vec event = as<vec>(Data["event"]);
    uvec idGK_fast = as<uvec>(Data["idGK_fast"]);
    mat W1 = as<mat>(Data["W1"]);
    mat W1s = as<mat>(Data["W1s"]);
    mat Wlong = as<mat>(Data["Wlong"]);
    mat Wlongs = as<mat>(Data["Wlongs"]);
    vec Pw = as<vec>(Data["Pw"]);
    // Initials
    vec Bs_gammas = as<vec>(initials["Bs_gammas"]);
    vec alphas = as<vec>(initials["alphas"]);
    vec phi_alphas = as<vec>(initials["phi_alphas"]);
    double tau_Bs_gammas = as<double>(initials["tau_Bs_gammas"]);
    double tau_alphas = as<double>(initials["tau_alphas"]);
    // Priors
    vec mean_Bs_gammas = as<vec>(priors["mean_Bs_gammas"]);
    mat Tau_Bs_gammas = as<mat>(priors["Tau_Bs_gammas"]);
    vec mean_alphas = as<vec>(priors["mean_alphas"]);
    mat Tau_alphas = as<mat>(priors["Tau_alphas"]);
    double A_tau_Bs_gammas = as<double>(priors["A_tau_Bs_gammas"]);
    double B_tau_Bs_gammas = as<double>(priors["B_tau_Bs_gammas"]);
    double Apost_tau_Bs_gammas = A_tau_Bs_gammas + 0.5 * as<double>(priors["rank_Tau_Bs_gammas"]);
    //
    bool shrink_alphas = as<bool>(priors["shrink_alphas"]);
    double A_tau_alphas = as<double>(priors["A_tau_alphas"]);
    double B_tau_alphas = as<double>(priors["B_tau_alphas"]);
    double Apost_tau_alphas = A_tau_alphas + 0.5 * as<double>(priors["rank_Tau_alphas"]);
    double A_phi_alphas = as<double>(priors["A_phi_alphas"]);
    double B_phi_alphas = as<double>(priors["B_phi_alphas"]);
    double Apost_phi_alphas = A_phi_alphas + 0.5;
    // Control parameters
    int n_iter = as<int>(control["n_iter"]);
    int n_burnin = as<int>(control["n_burnin"]);
    int n_block = as<int>(control["n_block"]);
    int n_thin = as<int>(control["n_thin"]);
    double target_acc = as<double>(control["target_acc"]);
    double eps1 = as<double>(control["eps1"]);
    double eps2 = as<double>(control["eps2"]);
    double eps3 = as<double>(control["eps3"]);
    double c0 = as<double>(control["c0"]);
    double mc1 = - as<double>(control["c1"]);
    bool adaptCov = as<bool>(control["adaptCov"]);
    // derived from control
    int total_it = n_iter + n_burnin;
    int its = ceil(total_it / n_block);
    int n_out = floor(n_iter / n_thin);
    int start_cov_update = floor(n_burnin / n_block) - 1;
    uvec keep = seqC(n_burnin + 1, total_it, n_thin);
    // initial
    int n_Bs_gammas = Bs_gammas.n_rows;
    int n_alphas = alphas.n_rows;
    double sigma_Bs_gammas = as<double>(scales["Bs_gammas"]);
    double sigma_alphas = as<double>(scales["alphas"]);
    mat Sigma_Bs_gammas = as<mat>(Covs["Bs_gammas"]);
    mat Sigma_alphas = as<mat>(Covs["alphas"]);
    // store results
    mat out_Bs_gammas(n_Bs_gammas, n_out, fill::zeros);
    mat out_alphas(n_alphas, n_out, fill::zeros);
    vec out_tau_Bs_gammas(n_out, fill::zeros);
    vec out_tau_alphas(n_out, fill::zeros);
    mat out_phi_alphas(n_alphas, n_out, fill::zeros);
    vec out_logPost(n_out, fill::zeros);
    RNGScope scope;
    int keep_iterator = 0;
    int check_iterator = 0;
    vec log_us = log(randu<vec>(its * n_block));
    double current_logpost = logPosterior_woRE_nogammas(event, W1, W1s, Bs_gammas,
                                               Wlong, Wlongs, alphas, Pw,
                                               mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                               mean_alphas, Tau_alphas, tau_alphas);
    for (int it = 0; it < its; ++it) {
        vec accept = vec(n_block, fill::zeros);
        mat block_Bs_gammas = mat(n_Bs_gammas, n_block, fill::zeros);
        mat block_alphas = mat(n_alphas, n_block, fill::zeros);
        mat rand_Bs_gammas = mvrnorm(n_block, sigma_Bs_gammas * Sigma_Bs_gammas);
        mat rand_alphas = mvrnorm(n_block, sigma_alphas * Sigma_alphas);
        vec new_Bs_gammas = vec(n_Bs_gammas, fill::zeros);
        vec new_alphas = vec(n_alphas, fill::zeros);
        for (int i = 0; i < n_block; ++i) {
            new_Bs_gammas = Bs_gammas + rand_Bs_gammas.col(i);
            new_alphas = alphas + rand_alphas.col(i);
            double new_logpost = logPosterior_woRE_nogammas(event, W1, W1s, new_Bs_gammas,
                                                   Wlong, Wlongs, new_alphas, Pw,
                                                   mean_Bs_gammas, Tau_Bs_gammas, tau_Bs_gammas,
                                                   mean_alphas, Tau_alphas, tau_alphas);
            double lRatio = new_logpost - current_logpost;
            int iter = it * n_block + i;
            if (log_us.at(iter) < lRatio) {
                accept.at(i) = 1;
                Bs_gammas = new_Bs_gammas;
                alphas = new_alphas;
                current_logpost = new_logpost;
            }
            block_Bs_gammas.col(i) = Bs_gammas;
            block_alphas.col(i) = alphas;
            // sample scales
            double Bpost_tau_Bs_gammas = B_tau_Bs_gammas + 0.5 * as_scalar(Bs_gammas.t() * Tau_Bs_gammas * Bs_gammas);
            tau_Bs_gammas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_Bs_gammas, 1.0 / Bpost_tau_Bs_gammas)));
            if (shrink_alphas) {
                for (int k2 = 0; k2 < n_alphas; ++k2) {
                    double Bpost_phi_alphas = B_phi_alphas + 0.5 * tau_alphas * pow(alphas.at(k2), 2);
                    phi_alphas.at(k2) = std::min(eps3,
                                  std::max(eps2, ::Rf_rgamma(Apost_phi_alphas, 1.0 / Bpost_phi_alphas)));
                }
                Tau_alphas.diag() = phi_alphas;
                double Bpost_tau_alphas = B_tau_alphas + 0.5 * as_scalar(alphas.t() * Tau_alphas * alphas);
                tau_alphas = std::min(eps3, std::max(eps2, ::Rf_rgamma(Apost_tau_alphas, 1.0 / Bpost_tau_alphas)));
            }
            // store results
            if ((unsigned)iter == keep[check_iterator]) {
                out_Bs_gammas.col(keep_iterator) = Bs_gammas;
                out_alphas.col(keep_iterator) = alphas;
                out_phi_alphas.col(keep_iterator) = phi_alphas;
                out_tau_Bs_gammas.at(keep_iterator) = tau_Bs_gammas;
                out_tau_alphas.at(keep_iterator) = tau_alphas;
                out_logPost.at(keep_iterator) = current_logpost;
                keep_iterator += 1;
                check_iterator += 1;
            }
        }
        double rate = mean(accept);
        double g1 = pow(it + 1, mc1);
        double g2 = c0 * g1;
        sigma_Bs_gammas = bounds_sigma(exp(log(sigma_Bs_gammas) + g2 * (rate - target_acc)), eps1, eps3);
        sigma_alphas = bounds_sigma(exp(log(sigma_alphas) + g2 * (rate - target_acc)), eps1, eps3);
        if (adaptCov && it > start_cov_update) {
            Sigma_Bs_gammas = bounds_Cov(Sigma_Bs_gammas + g1 * (cov(block_Bs_gammas.t()) - Sigma_Bs_gammas), eps2, eps3);
            Sigma_alphas = bounds_Cov(Sigma_alphas + g1 * (cov(block_alphas.t()) - Sigma_alphas), eps2, eps3);
        }
    }
    // export
    return List::create(
        Named("mcmc") = List::create(
            Named("Bs_gammas") = out_Bs_gammas,
            Named("alphas") = out_alphas,
            Named("phi_alphas") = out_phi_alphas,
            Named("tau_Bs_gammas") = out_tau_Bs_gammas,
            Named("tau_alphas") = out_tau_alphas
        ),
        
        Named("logPost") = out_logPost,
        
        Named("scales") = List::create(
            Named("sigma") = List::create(
                Named("Bs_gammas") = sigma_Bs_gammas,
                Named("alphas") = sigma_alphas
            ),
            Named("Sigma") = List::create(
                Named("Bs_gammas") = Sigma_Bs_gammas,
                Named("alphas") = Sigma_alphas
            )
        )
    );
}
