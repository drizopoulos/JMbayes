#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double logPosterior (double temp, arma::vec event, arma::uvec idGK_fast,
                     arma::mat W1, arma::mat W1s, arma::vec Bs_gammas,
                     arma::mat W2, arma::mat W2s, arma::vec gammas,
                     arma::mat Wlong, arma::mat Wlongs, arma::vec alphas,
                     arma::vec Pw,
                     arma::vec mean_Bs_gammas, arma::mat Tau_Bs_gammas,
                     arma::vec mean_gammas, arma::mat Tau_gammas,
                     arma::vec mean_alphas, arma::mat Tau_alphas,
                     double tauBs, double A_tauBs, double B_tauBs) {
    // log-likelihood
    vec log_h = W1 * Bs_gammas + W2 * gammas + Wlong * alphas;
    vec H = Pw % exp(W1s * Bs_gammas + W2s * gammas + Wlongs * alphas);
    //
    double out = temp * (sum(event % log_h) - sum(H)) -
        0.5 * tauBs * as_scalar((Bs_gammas - mean_Bs_gammas).t() * Tau_Bs_gammas * (Bs_gammas - mean_Bs_gammas)) -
        0.5 * as_scalar((gammas - mean_gammas).t() * Tau_gammas * (gammas - mean_gammas)) -
        0.5 * as_scalar((alphas - mean_alphas).t() * Tau_alphas * (alphas - mean_alphas)) +
        (A_tauBs - 1) * log(tauBs) - B_tauBs * tauBs;
    return(out);
}

// [[Rcpp::export]]
double logPosterior_nogammas (double temp, arma::vec event, arma::uvec idGK_fast,
                     arma::mat W1, arma::mat W1s, arma::vec Bs_gammas,
                     arma::mat Wlong, arma::mat Wlongs, arma::vec alphas,
                     arma::vec Pw,
                     arma::vec mean_Bs_gammas, arma::mat Tau_Bs_gammas,
                     arma::vec mean_alphas, arma::mat Tau_alphas,
                     double tauBs, double A_tauBs, double B_tauBs) {
    // log-likelihood
    vec log_h = W1 * Bs_gammas + Wlong * alphas;
    vec H = Pw % exp(W1s * Bs_gammas + Wlongs * alphas);
    //
    double out = temp * (sum(event % log_h) - sum(H)) -
        0.5 * tauBs * as_scalar((Bs_gammas - mean_Bs_gammas).t() * Tau_Bs_gammas * (Bs_gammas - mean_Bs_gammas)) -
        0.5 * as_scalar((alphas - mean_alphas).t() * Tau_alphas * (alphas - mean_alphas)) +
        (A_tauBs - 1) * log(tauBs) - B_tauBs * tauBs;
    return(out);
}

