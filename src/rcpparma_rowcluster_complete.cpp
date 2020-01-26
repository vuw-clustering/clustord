// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// [[Rcpp::export]]
double rcpparma_Rclusterll(arma::mat & ymat, arma::cube & theta,
                           arma::mat & pprm, arma::colvec piv,
                           int RG, int p, int n, bool partial) {
    double llc = 0;
    int ymatiijj;

    for (int rr=0; rr<RG; ++rr) {
        arma::mat log_thetaymat(n,p);
        log_thetaymat.fill(0);

        for (int jj=0; jj<p; ++jj) {
            for (int ii=0; ii<n; ++ii) {
                if (arma::is_finite(ymat(ii,jj))) {
                    ymatiijj = ymat(ii,jj)-1;
                    log_thetaymat(ii,jj)=log(theta(rr,jj,ymatiijj));
                }
            }
        }
        llc += sum(pprm.col(rr).t()*log_thetaymat);
    }

    if (!partial) llc += sum(pprm*log(piv));

    return llc;
}