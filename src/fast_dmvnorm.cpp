#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec dmvnorm2(const arma::mat& x, const arma::rowvec& mean, 
                   const arma::mat& sigma, bool logd = false) {
    int n = x.n_rows;
    int xdim = x.n_cols;
    arma::vec out(n);
    arma::mat rooti = arma::trans(arma::inv(trimatu(arma::chol(sigma))));
    double rootisum = arma::sum(log(rooti.diag()));
    double log2pi = std::log(2.0 * M_PI);
    double constants = -(static_cast<double>(xdim)/2.0) * log2pi;
    for (int i = 0; i < n; ++i) {
        arma::vec z = rooti * arma::trans(x.row(i) - mean) ;
        out(i) = constants - 0.5 * arma::sum(z%z) + rootisum;
    }
    if (logd == false) {
        out = exp(out);
    }
    return(out);
}

