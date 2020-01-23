// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
double rcpparma_updateliC(double li, const arma::colvec & columnclusters,
                          const int RG, const int n, const int p,
                          const arma::mat & ymat, const arma::cube & theta,
                          const arma::colvec & pivec, const arma::rowvec & kappavec) {

    double alpha = 1;
    int cc;
    arma::mat Aair(n,RG);
    arma::colvec maxAair(n);
    double logDaim;
    double Eav;

    // Initialize arrays to NA
    Aair.fill(NA_REAL);
    maxAair.fill(NA_REAL);

    for(int jj=0; jj<p; ++jj) {
        cc = columnclusters[jj]-1;
        alpha = alpha*kappavec[cc];
    }
    // printf("alpha = %f\n",alpha);
    if(alpha > 0) {
        double pirr;
        double ymatiijj;
        double thetarrcciijj;

        for (int ii=0; ii<n; ++ii) {
            for (int rr=0; rr<RG; ++rr) {
                pirr = pivec[rr];
                Aair(ii,rr) = log(pirr);
                for (int jj=0; jj<p; ++jj) {
                    cc = columnclusters[jj]-1;
                    ymatiijj=ymat(ii,jj)-1;
                    thetarrcciijj = theta(rr,cc,ymatiijj);
                    if (thetarrcciijj > 0) {
                        Aair(ii,rr) += log(thetarrcciijj);
                        // printf("Aair = %f\n",Aair(ii,rr));
                    }
                }
            }
        }
    }

    maxAair = max(Aair,1);
    for (int ii=0; ii<n; ++ii) {
        // printf("maxAair = %f\n",maxAair(ii));
        logDaim = log(sum(exp(Aair.row(ii)-maxAair(ii))));
        if (arma::is_finite(logDaim)) Eav += maxAair(ii) + logDaim;
    }

    Eav += log(alpha);
    // printf("Eav = %f\n",Eav);
    if (arma::is_finite(Eav)) {
        li += exp(Eav);
    }

    // arma::colvec a(2);
    // a(0) = NA_REAL;
    // printf("max of NA and 0 is %f\n",max(a));
    // printf("sum of NA and 0 is %f\n",sum(a));

    return li;
}

// [[Rcpp::export]]
double rcpparma_updateliR(double li, const arma::rowvec & rowclusters,
                          const int CG, const int n, const int p,
                          const arma::mat & ymat, const arma::cube & theta,
                          const arma::colvec & pivec, const arma::rowvec & kappavec) {

    double alpha = 1;
    int rr;
    arma::mat Aajc(p,CG);
    arma::colvec maxAajc(p);
    double logDajm;
    double Eav;

    // Initialize arrays to NA
    Aajc.fill(NA_REAL);
    maxAajc.fill(NA_REAL);

    for(int ii=0; ii<n; ++ii) {
        rr = rowclusters[ii]-1;
        alpha = alpha*pivec[rr];
    }
    // printf("alpha = %f\n",alpha);
    if(alpha > 0) {
        double kappacc;
        double ymatiijj;
        double thetarrcciijj;

        for (int jj=0; jj<p; ++jj) {
            for (int cc=0; cc<CG; ++cc) {
                kappacc = kappavec[cc];
                Aajc(jj,cc) = log(kappacc);
                for (int ii=0; ii<n; ++ii) {
                    rr = rowclusters[ii]-1;
                    ymatiijj=ymat(ii,jj)-1;
                    thetarrcciijj = theta(rr,cc,ymatiijj);
                    if (thetarrcciijj > 0) {
                        Aajc(jj,cc) += log(thetarrcciijj);
                        // printf("Aajc = %f\n",Aajc(jj,cc));
                    }
                }
            }
        }
    }

    maxAajc = max(Aajc,1);
    for (int jj=0; jj<p; ++jj) {
        // printf("maxAajc = %f\n",maxAajc(jj));
        logDajm = log(sum(exp(Aajc.row(jj)-maxAajc(jj))));
        if (arma::is_finite(logDajm)) Eav += maxAajc(jj) + logDajm;
    }

    Eav += log(alpha);
    // printf("Eav = %f\n",Eav);
    if (arma::is_finite(Eav)) {
        li += exp(Eav);
    }

    // arma::colvec a(2);
    // a(0) = NA_REAL;
    // printf("max of NA and 0 is %f\n",max(a));
    // printf("sum of NA and 0 is %f\n",sum(a));

    return li;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
rcpparma_updateliC(0,c(1,1,1),2,2,3,matrix(c(1,1,2,2,1,1),nrow=2),
                    array(1:12,dim=c(2,3,2)),c(0.5,0.5),c(0.2,0.8))
*/
