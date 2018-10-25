#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

using namespace Rcpp;
using namespace arma;

arma::vec rowsum_svft(const vec& x_, const uvec& group) {
    vec x = cumsum(x_);
    vec out = x.elem(group);
    out.insert_rows(0, 1);
    unsigned int n = out.n_elem;
    out = out.rows(1, n - 1) - out.rows(0, n - 2);
    return(out);
}

arma::vec Vpnorm_svft(const vec& x) {
    int n = x.size();
    vec res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = ::Rf_pnorm5(x[i], 0.0, 1.0, 1, 0);
    }
    return(res);
}

arma::field<arma::vec> List2Field_vec_svft(const List & Vecs) {
    int n_list = Vecs.size();
    arma::field<arma::vec> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<arma::vec>(Vecs[i]);
    }
    return(res);
}

arma::field<arma::uvec> List2Field_uvec_svft(const List & uVecs, bool substract1 = true) {
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

arma::field<arma::mat> List2Field_mat_svft(const List & Mats) {
    int n_list = Mats.size();
    arma::field<arma::mat> res(n_list);
    for (int i = 0; i < n_list; ++i) {
        res.at(i) = as<arma::mat>(Mats[i]);
    }
    return(res);
}

arma::field<arma::vec> lin_predF_svft(const field<vec>& Xbetas, const field<mat>& Z,
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

arma::mat lin_pred_matF_svft(const field<vec>& Xbetas, const field<mat>& Z,
                             const mat& b, const field<mat>& U,
                             const field<uvec>& RE_inds, const field<uvec>& id,
                             const field<uvec>& col_inds, const uvec& row_inds,
                             const int& nrows, const int& ncols,
                             const CharacterVector& trans_Funs) {
    int n_field = Xbetas.size();
    mat out = mat(nrows, ncols, fill::zeros);
    for (int i = 0; i < n_field; ++i) {
        mat bb = b.cols(RE_inds.at(i));
        vec Zb = sum(Z.at(i) % bb.rows(id.at(i)), 1);
        if (trans_Funs[i] == "identity") {
            out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % (Xbetas.at(i) + Zb);
        } else if (trans_Funs[i] == "expit") {
            vec exp_eta = exp(Xbetas.at(i) + Zb);
            out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % (exp_eta / (1 + exp_eta));
        } else if (trans_Funs[i] == "exp") {
            out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % exp(Xbetas.at(i) + Zb);
        } else if (trans_Funs[i] == "log") {
            out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % log(Xbetas.at(i) + Zb);
        } else if (trans_Funs[i] == "log2") {
            out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % log2(Xbetas.at(i) + Zb);
        } else if (trans_Funs[i] == "log10") {
            out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % log10(Xbetas.at(i) + Zb);
        } else if (trans_Funs[i] == "sqrt") {
            out.submat(row_inds, col_inds.at(i)) = U.at(i).each_col() % sqrt(Xbetas.at(i) + Zb);
        }
    }
    return(out);
}

arma::vec log_longF_svft(const field<vec>& y, const field<vec>& eta,
                         const CharacterVector& fams, const CharacterVector& links,
                         const List& sigmas, const field<uvec>& id, const int& n) {
    int n_outcomes = y.size();
    vec out(n, fill::zeros);
    for (int i = 0; i < n_outcomes; ++i) {
        uvec id_i = id.at(i);
        vec y_i = y.at(i);
        vec eta_i = eta.at(i);
        if (y_i.size() == 0) {
            out += 0;
        } else {
            if (fams[i] == "gaussian") {
                double sigma_i = as<double>(sigmas[i]);
                vec log_dens = - 0.5 * pow((y_i - eta_i) / sigma_i, 2);
                out += rowsum_svft(log_dens, id_i);
            } else if (fams[i] == "binomial") {
                if (links[i] == "logit") {
                    vec pr = exp(eta_i) / (1 + exp(eta_i));
                    vec log_dens = y_i % log(pr) + (1 - y_i) % log(1 - pr);
                    out += rowsum_svft(log_dens, id_i);
                } else if (links[i] == "probit") {
                    vec pr = Vpnorm_svft(eta_i);
                    vec log_dens = y_i % log(pr) + (1 - y_i) % log(1 - pr);
                    out += rowsum_svft(log_dens, id_i);
                } else if (links[i] == "cloglog") {
                    vec pr = - exp(- exp(eta_i)) + 1;
                    vec log_dens = y_i % log(pr) + (1 - y_i) % log(1 - pr);
                    out += rowsum_svft(log_dens, id_i);
                }
            } else if (fams[i] == "poisson") {
                vec mu = exp(eta_i);
                vec log_dens = y_i % log(mu) - mu;
                out += rowsum_svft(log_dens, id_i);
            }
        }
    }
    return(out);
}

arma::vec log_postREF_svft(const mat& b_mat, const vec& Bs_gammas, const vec& gammas,
                           const vec& alphas, const field<vec>& y, const field<vec>& Xbetas,
                           const field<mat>& Z,
                           const field<uvec>& RE_inds, const field<uvec>& RE_inds2,
                           const field<uvec>& idL, const field<uvec>& idL2,
                           const CharacterVector& fams, const CharacterVector& links,
                           const List& sigmas, const mat& invD,
                           const int& n, const int& ns, const int& n_alphas,
                           const mat& W1s,
                           const mat& W2s,
                           const field<vec>& XXsbetas,
                           const field<mat>& ZZs,
                           const field<mat>& Us, const vec& Pw, const uvec& idGK,
                           const field<uvec>& idTs,
                           const field<uvec>& col_inds,
                           const uvec& row_inds_Us,
                           const CharacterVector& trans_Funs) {
    field<vec> eta = lin_predF_svft(Xbetas, Z, b_mat, RE_inds, idL);
    vec log_pyb = log_longF_svft(y, eta, fams, links, sigmas, idL2, n);
    vec log_pb = - 0.5 * sum((b_mat * invD) % b_mat, 1);
    mat Wlongs = lin_pred_matF_svft(XXsbetas, ZZs, b_mat, Us, RE_inds2, idTs,
                                    col_inds, row_inds_Us, ns, n_alphas, trans_Funs);
    vec H = rowsum_svft(Pw % exp(W1s * Bs_gammas + W2s * gammas + Wlongs * alphas), idGK);
    vec log_ptb = - H;
    vec out = log_pyb + log_ptb + log_pb;
    return(out);
}

arma::vec survPred_svft(const mat& b_mat, const vec& Bs_gammas, const vec& gammas,
                        const vec& alphas,
                        const field<uvec>& RE_inds2,
                        const int& ns, const int& n_alphas,
                        const mat& W1s,
                        const mat& W2s,
                        const field<vec>& XXsbetas,
                        const field<mat>& ZZs,
                        const field<mat>& Us, const vec& Pw, const uvec& idGK,
                        const field<uvec>& idTs,
                        const field<uvec>& col_inds,
                        const uvec& row_inds_Us,
                        const CharacterVector& trans_Funs) {
    mat Wlongs = lin_pred_matF_svft(XXsbetas, ZZs, b_mat, Us, RE_inds2, idTs,
                                    col_inds, row_inds_Us, ns, n_alphas, trans_Funs);
    vec H = rowsum_svft(Pw % exp(W1s * Bs_gammas + W2s * gammas + Wlongs * alphas), idGK);
    vec log_ptb = - H;
    vec out = log_ptb;
    return(out);
}

// [[Rcpp::export]]
arma::vec log_post_RE_svft(arma::vec b, List Data) {
    mat b_mat = conv_to<rowvec>::from(b);
    vec Bs_gammas = as<vec>(Data["Bs_gammas"]);
    vec gammas = as<vec>(Data["gammas"]);
    vec alphas = as<vec>(Data["alphas"]);
    List y = as<List>(Data["y"]);
    field<vec> yF = List2Field_vec_svft(y);
    List Xbetas = as<List>(Data["Xbetas"]);
    field<vec> XbetasF = List2Field_vec_svft(Xbetas);
    List Z = as<List>(Data["Z"]);
    field<mat> ZF = List2Field_mat_svft(Z);
    List RE_inds = as<List>(Data["RE_inds"]);
    field<uvec> RE_indsF = List2Field_uvec_svft(RE_inds);
    List RE_inds2 = as<List>(Data["RE_inds2"]);
    field<uvec> RE_inds2F = List2Field_uvec_svft(RE_inds2);
    List idL = as<List>(Data["idL"]);
    field<uvec> idLF = List2Field_uvec_svft(idL);
    List idL2 = as<List>(Data["idL2"]);
    field<uvec> idL2F = List2Field_uvec_svft(idL2);
    CharacterVector fams = as<CharacterVector>(Data["fams"]);
    CharacterVector links = as<CharacterVector>(Data["links"]);
    List sigmas = as<List>(Data["sigmas"]);
    mat invD = as<mat>(Data["invD"]);
    int n = b_mat.n_rows;
    int n_alphas = alphas.n_rows;
    mat W1s = as<mat>(Data["W1s"]);
    mat W2s = as<mat>(Data["W2s"]);
    List XXsbetas = as<List>(Data["XXsbetas"]);
    field<vec> XXsbetasF = List2Field_vec_svft(XXsbetas);
    List ZZs = as<List>(Data["ZZs"]);
    field<mat> ZZsF = List2Field_mat_svft(ZZs);
    List Us = as<List>(Data["Us"]);
    field<mat> UsF = List2Field_mat_svft(Us);
    vec Pw = as<vec>(Data["Pw"]);
    int ns = Pw.n_rows;
    uvec idGK = as<uvec>(Data["idGK"]);
    List idTs = as<List>(Data["idTs"]);
    field<uvec> idTsF = List2Field_uvec_svft(idTs);
    List col_inds = as<List>(Data["col_inds"]);
    field<uvec> col_indsF = List2Field_uvec_svft(col_inds);
    uvec row_inds_Us = as<uvec>(Data["row_inds_Us"]) - 1;
    CharacterVector trans_Funs = as<CharacterVector>(Data["trans_Funs"]);
    arma::vec out = log_postREF_svft(b_mat, Bs_gammas, gammas, alphas, yF, XbetasF, 
                                     ZF, RE_indsF, RE_inds2F, idLF, idL2F, fams, links, 
                                     sigmas, invD, n, ns, n_alphas, W1s, W2s, XXsbetasF, 
                                     ZZsF, UsF, Pw, idGK, idTsF, col_indsF, row_inds_Us, 
                                     trans_Funs);
    return(out);
}

// [[Rcpp::export]]
arma::vec survPred_svft_2(arma::vec b, List Data) {
    mat b_mat = conv_to<rowvec>::from(b);
    vec Bs_gammas = as<vec>(Data["Bs_gammas"]);
    vec gammas = as<vec>(Data["gammas"]);
    vec alphas = as<vec>(Data["alphas"]);
    List RE_inds2 = as<List>(Data["RE_inds2"]);
    field<uvec> RE_inds2F = List2Field_uvec_svft(RE_inds2);
    mat W1s = as<mat>(Data["W1s"]);
    mat W2s = as<mat>(Data["W2s"]);
    List XXsbetas = as<List>(Data["XXsbetas"]);
    field<vec> XXsbetasF = List2Field_vec_svft(XXsbetas);
    List ZZs = as<List>(Data["ZZs"]);
    field<mat> ZZsF = List2Field_mat_svft(ZZs);
    List Us = as<List>(Data["Us"]);
    field<mat> UsF = List2Field_mat_svft(Us);
    vec Pw = as<vec>(Data["Pw"]);
    int ns = Pw.n_rows;
    int n_alphas = alphas.n_rows;
    uvec idGK = as<uvec>(Data["idGK"]);
    List idTs = as<List>(Data["idTs"]);
    field<uvec> idTsF = List2Field_uvec_svft(idTs);
    List col_inds = as<List>(Data["col_inds"]);
    field<uvec> col_indsF = List2Field_uvec_svft(col_inds);
    uvec row_inds_Us = as<uvec>(Data["row_inds_Us"]) - 1;
    CharacterVector trans_Funs = as<CharacterVector>(Data["trans_Funs"]);
    arma::vec out = survPred_svft(b_mat, Bs_gammas, gammas, alphas, RE_inds2F, 
                                  ns, n_alphas, W1s, W2s, XXsbetasF, 
                                  ZZsF, UsF, Pw, idGK, idTsF, col_indsF, row_inds_Us, 
                                  trans_Funs);
    return(out);
}
