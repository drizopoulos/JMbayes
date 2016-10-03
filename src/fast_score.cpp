#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec gradient_logPosterior (double temp, arma::vec event, arma::uvec idGK_fast,
                           arma::vec event_colSumsW1, arma::mat W1s, arma::vec Bs_gammas,
                           arma::vec event_colSumsW2, arma::mat W2s, arma::vec gammas,
                           arma::vec event_colSumsWlong, arma::mat Wlongs, arma::vec alphas,
                           arma::vec Pw,
                           arma::vec mean_Bs_gammas, arma::mat Tau_Bs_gammas,
                           arma::vec mean_gammas, arma::mat Tau_gammas,
                           arma::vec mean_alphas, arma::mat Tau_alphas,
                           double tauBs, double A_tauBs, double B_tauBs) {
    vec Int = Pw % exp(W1s * Bs_gammas + W2s * gammas + Wlongs * alphas);
    // gradient Bs_gammas
    unsigned int n_Bs_gammas = Bs_gammas.n_rows;
    vec score_Bs_gammas = vec(n_Bs_gammas, fill::zeros);
    for (uword i = 0; i < n_Bs_gammas; ++i) {
        vec H = W1s.col(i) % Int;
        score_Bs_gammas[i] = temp * (event_colSumsW1[i] - sum(H));
    }
    score_Bs_gammas = score_Bs_gammas - tauBs * (Tau_Bs_gammas * (Bs_gammas - mean_Bs_gammas));
    // gradient gammas
    unsigned int n_gammas = gammas.n_rows;
    vec score_gammas = vec(n_gammas, fill::zeros);
    for (uword i = 0; i < n_gammas; ++i) {
        vec H = W2s.col(i) % Int;
        score_gammas[i] = temp * (event_colSumsW2[i] - sum(H));
    }
    score_gammas = score_gammas - (Tau_gammas * (gammas - mean_gammas));
    // gradient alphas
    unsigned int n_alphas = alphas.n_rows;
    vec score_alphas = vec(n_alphas, fill::zeros);
    for (uword i = 0; i < n_alphas; ++i) {
        vec H = Wlongs.col(i) % Int;
        score_alphas[i] = temp * (event_colSumsWlong[i] - sum(H));
    }
    score_alphas = score_alphas - (Tau_alphas * (alphas - mean_alphas));
    // gradient tauBs
    vec z = Bs_gammas - mean_Bs_gammas;
    vec score_tauBs = tauBs * (- 0.5 * (z.t() * Tau_Bs_gammas * z) + (A_tauBs - 1) / tauBs - B_tauBs);
    // export results
    vec out = join_cols(score_Bs_gammas, score_gammas);
    out = join_cols(out, score_alphas);
    out = join_cols(out, score_tauBs);
    return(out);
}

// [[Rcpp::export]]
arma::vec gradient_logPosterior_nogammas (double temp, arma::vec event, arma::uvec idGK_fast,
                                 arma::vec event_colSumsW1, arma::mat W1s, arma::vec Bs_gammas,
                                 arma::vec event_colSumsWlong, arma::mat Wlongs, arma::vec alphas,
                                 arma::vec Pw,
                                 arma::vec mean_Bs_gammas, arma::mat Tau_Bs_gammas,
                                 arma::vec mean_alphas, arma::mat Tau_alphas,
                                 double tauBs, double A_tauBs, double B_tauBs) {
    vec Int = Pw % exp(W1s * Bs_gammas + Wlongs * alphas);
    // gradient Bs_gammas
    unsigned int n_Bs_gammas = Bs_gammas.n_rows;
    vec score_Bs_gammas = vec(n_Bs_gammas, fill::zeros);
    for (uword i = 0; i < n_Bs_gammas; ++i) {
        vec H = W1s.col(i) % Int;
        score_Bs_gammas[i] = temp * (event_colSumsW1[i] - sum(H));
    }
    score_Bs_gammas = score_Bs_gammas - tauBs * (Tau_Bs_gammas * (Bs_gammas - mean_Bs_gammas));
    // gradient alphas
    unsigned int n_alphas = alphas.n_rows;
    vec score_alphas = vec(n_alphas, fill::zeros);
    for (uword i = 0; i < n_alphas; ++i) {
        vec H = Wlongs.col(i) % Int;
        score_alphas[i] = temp * (event_colSumsWlong[i] - sum(H));
    }
    score_alphas = score_alphas - (Tau_alphas * (alphas - mean_alphas));
    // gradient tauBs
    vec z = Bs_gammas - mean_Bs_gammas;
    vec score_tauBs = tauBs * (- 0.5 * (z.t() * Tau_Bs_gammas * z) + (A_tauBs - 1) / tauBs - B_tauBs);
    // export results
    vec out = join_cols(score_Bs_gammas, score_alphas);
    out = join_cols(out, score_tauBs);
    return(out);
}
