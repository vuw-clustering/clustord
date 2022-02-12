// -*- mode: C++; cc-indent-level: 4; cc-basic-offset: 4; indent-tabs-mode: nil; -*-

#include "Rcpp.h"

using namespace Rcpp;

double rcpp_expit(double x) {

    double out = 1/(1 + exp(-x));
    return out;
}

double rcpp_logit(double x) {

    double out = log(x/(1-x));
    return out;
}

void rcpp_unpack(const String & model,
                 const NumericVector & invect,
                 const IntegerVector & paramlengths,
                 NumericVector & mu,
                 NumericVector & phi,
                 NumericVector & rowc_coef,
                 NumericVector & colc_coef,
                 NumericMatrix & rowc_colc_coef,
                 NumericVector & row_coef,
                 NumericVector & col_coef,
                 NumericMatrix & rowc_col_coef,
                 NumericMatrix & colc_row_coef,
                 NumericMatrix & rowc_cov_coef,
                 NumericMatrix & colc_cov_coef,
                 NumericVector & cov_coef,
                 const int & RG, const int & CG,
                 const int & p, const int & n, const int & q,
                 const bool constraint_sum_zero) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: paramlengths["rowc"], else you will get a
    // fatal error that will crash R;

    int ind, kk, rr, cc, i, j, l;
    int nelts = 0;

    if (model == "OSM") {
        mu[0] = 0;
        for (kk=1; kk < q; kk++) {
            ind = kk-1;
            // Rcout << "The index of invect : " << ind << "\n";
            mu[kk] = invect[ind];
        }
        nelts += q-1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of mu : " << mu << "\n";

        // Convert to phi from u, where u can vary between -Inf and +Inf but phi
        // must be between 0 and 1, and phi_k >= phi_k-1
        NumericVector u (q-1);
        for (kk=1; kk < q-1; kk++) {
            ind = nelts + kk-1;
            // Rcout << "The index of invect : " << ind << "\n";
            u[kk] = invect[ind];
        }
        // Rcout << "The value of u : " << u << "\n";

        phi[0] = 0;
        phi[1] = rcpp_expit(u[1]);
        if (q == 3) {
            phi[2] = 1;
        } else {
            double phiraw = 0;
            for (kk=2; kk < q-1; kk++) {
                phiraw += exp(u[kk]);
                // Rcout << "The value of phiraw : " << phiraw << "\n";
                phi[kk] = rcpp_expit(u[2] + phiraw);
            }
        }
        phi[q-1] = 1;
        nelts += q-2;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of phi : " << phi << "\n";

    } else if (model == "POM") {
        for (kk=0; kk < q-1; kk++) {
            mu[kk] = invect[kk];
        }
        mu.sort();
        nelts += q-1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of mu : " << mu << "\n";
    } else if (model == "Binary") {
        mu[0] = invect[0];
        nelts += 1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of mu : " << mu << "\n";
    }

    if (paramlengths["rowc"] > 0) {
        if (constraint_sum_zero) {
            rowc_coef[RG-1] = 0;
            for (rr=0; rr < RG-1; rr++) {
                ind = nelts + rr;
                // Rcout << "The index of invect : " << ind << "\n";
                rowc_coef[rr] = invect[ind];
                rowc_coef[RG-1] -= invect[ind];
            }
        } else {
            rowc_coef[0] = 0;
            for (rr=1; rr < RG; rr++) {
                ind = nelts + rr-1;
                // Rcout << "The index of invect : " << ind << "\n";
                rowc_coef[rr] = invect[ind];
            }
        }
        nelts += RG-1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of rowc_coef : " << rowc_coef << "\n";
    }
    if (paramlengths["colc"] > 0) {
        if (constraint_sum_zero) {
            colc_coef[CG-1] = 0;
            for (cc=0; cc < CG-1; cc++) {
                ind = nelts + cc;
                // Rcout << "The index of invect : " << ind << "\n";
                colc_coef[cc] = invect[ind];
                colc_coef[CG-1] -= invect[ind];
            }
        } else {
            colc_coef[0] = 0;
            for (cc=1; cc < CG; cc++) {
                ind = nelts + cc-1;
                // Rcout << "The index of invect : " << ind << "\n";
                colc_coef[cc] = invect[ind];
            }
        }
        nelts += CG-1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of colc_coef : " << colc_coef << "\n";
    }
    if (paramlengths["rowc_colc"] > 0) {
        for (cc=0; cc < CG; cc++) {
            rowc_colc_coef(RG-1,cc) = 0;
        }
        for (rr=0; rr < RG-1; rr++) {
            rowc_colc_coef(rr,0) = 0;
            for (cc=1; cc < CG; cc++) {
                ind = nelts + rr*(CG-1)+cc-1;
                // Rcout << "The index of invect : " << ind << "\n";
                rowc_colc_coef(rr,0) -= invect[ind];
                rowc_colc_coef(RG-1,0) += invect[ind];
                rowc_colc_coef(rr,cc) = invect[ind];
                rowc_colc_coef(RG-1,cc) -= invect[ind];
            }
        }
        nelts += (RG-1)*(CG-1);
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of rowc_colc_coef : " << rowc_colc_coef << "\n";
    }

    if (paramlengths["row"] > 0) {
        for (i=0; i < n; i++) {
            ind = nelts + i;
            // Rcout << "The index of invect : " << ind << "\n";
            row_coef[i] = invect[ind];
        }
        nelts += n;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of row_coef : " << row_coef << "\n";
    }
    if (paramlengths["col"] > 0) {
        for (j=0; j < p; j++) {
            ind = nelts + j;
            // Rcout << "The index of invect : " << ind << "\n";
            col_coef[j] = invect[ind];
        }
        nelts += p;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of col_coef : " << col_coef << "\n";
    }

    if (paramlengths["rowc_col"] > 0) {
        // NOTE: have to make sure to fill the matrix of rowc_col coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        for (rr=0; rr < RG; rr++) {
            for (j=0; j < p; j++) {
                rowc_col_coef(rr,j) = invect[nelts + rr*p + j];
            }
        }
        nelts += RG*p;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of rowc_col_coef : " << rowc_col_coef << "\n";
    }
    if (paramlengths["colc_row"] > 0) {
        // NOTE: have to make sure to fill the matrix of colc_row coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        for (cc=0; cc < CG; cc++) {
            for (i=0; i < n; i++) {
                colc_row_coef(cc,i) = invect[nelts + cc*n + i];
            }
        }
        nelts += CG*n;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of colc_row_coef : " << colc_row_coef << "\n";
    }

    if (paramlengths["rowc_cov"] > 0 && RG > 0) {
        // NOTE: have to make sure to fill the matrix of rowc_cov coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        int nrowccov = paramlengths["rowc_cov"]/RG;
        // Rcout << "The value of nrowccov : " << nrowccov << "\n";
        for (rr=0; rr < RG; rr++) {
            for (l=0; l < nrowccov; l++) {
                ind = nelts + rr*nrowccov + l;
                // Rcout << "The index of invect : " << ind << "\n";
                rowc_cov_coef(rr,l) = invect[ind];
            }
        }
        nelts += RG*nrowccov;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of rowc_cov_coef : " << rowc_cov_coef << "\n";
    }
    if (paramlengths["colc_cov"] > 0 && CG > 0) {
        // NOTE: have to make sure to fill the matrix of colc_cov coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        int ncolccov = paramlengths["colc_cov"]/CG;
        for (cc=0; cc < CG; cc++) {
            for (l=0; l < ncolccov; l++) {
                ind = nelts + cc*ncolccov + l;
                // Rcout << "The index of invect : " << ind << "\n";
                colc_cov_coef(cc,l) = invect[ind];
            }
        }
        nelts += CG*ncolccov;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of colc_cov_coef : " << colc_cov_coef << "\n";
    }

    int ncov = paramlengths["cov"];
    if (ncov > 0) {
        for (l=0; l < ncov; l++) {
            ind = nelts + l;
            // Rcout << "The index of invect : " << ind << "\n";
            cov_coef[l] = invect[ind];
        }
        nelts += ncov;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of cov_coef : " << cov_coef << "\n";
    }
}

double rcpp_linear_part(const NumericMatrix ydf,
                        const NumericMatrix & rowcmm,
                        const NumericMatrix & colcmm,
                        const NumericMatrix & covmm,
                        const IntegerVector & paramlengths,
                        const NumericVector & rowc_coef,
                        const NumericVector & colc_coef,
                        const NumericMatrix & rowc_colc_coef,
                        const NumericVector & row_coef,
                        const NumericVector & col_coef,
                        const NumericMatrix & rowc_col_coef,
                        const NumericMatrix & colc_row_coef,
                        const NumericMatrix & rowc_cov_coef,
                        const NumericMatrix & colc_cov_coef,
                        const NumericVector & cov_coef,
                        const int & RG, const int & CG,
                        const int & p, const int & n, const int & q,
                        const int & nrowccov, const int & ncolccov,
                        const int & rr, const int & cc,
                        const int & ij, const int & ii,
                        const int & jj) {

    double linear_part = 0;
    int ll;

    // Rcout << "The R-based value of ii : " << ii+1 << "\n";
    // Rcout << "The R-based value of jj : " << jj+1 << "\n";

    if (paramlengths["rowc"] > 0) {
        linear_part += rowc_coef[rr];
        // Rcout << "The value of linear_part with rowc: " << linear_part << "\n";
    }
    if (paramlengths["colc"] > 0) {
        linear_part += colc_coef[cc];
        // Rcout << "The value of linear_part with colc: " << linear_part << "\n";
    }
    if (paramlengths["rowc_colc"] > 0) {
        linear_part += rowc_colc_coef(rr,cc);
        // Rcout << "The value of linear_part with rowc_colc: " << linear_part << "\n";
    }

    if (paramlengths["row"] > 0) {
        linear_part += row_coef[ii];
        // Rcout << "The value of linear_part with row: " << linear_part << "\n";
    }
    if (paramlengths["col"] > 0) {
        linear_part += col_coef[jj];
        Rcout << "The value of col_coef: " << col_coef[jj] << "\n";
        // Rcout << "The value of linear_part with col: " << linear_part << "\n";
    }
    if (paramlengths["rowc_col"] > 0) {
        linear_part += rowc_col_coef(rr,jj);
        // Rcout << "The value of linear_part with rowc_col: " << linear_part << "\n";
    }
    if (paramlengths["colc_row"] > 0) {
        linear_part += colc_row_coef(cc,ii);
        // Rcout << "The value of linear_part with colc_row: " << linear_part << "\n";
    }

    if (paramlengths["rowc_cov"] > 0) {
        // NumericVector temp = rowcmm.row(ij);
        // Rcout << "The value of rowcmm: " << temp << "\n";

        for (ll=0; ll < nrowccov; ll++) {
            linear_part += rowcmm(ij,ll)*rowc_cov_coef(rr,ll);
        }
        // Rcout << "The value of linear_part with rowc_cov: " << linear_part << "\n";
    }
    if (paramlengths["colc_cov"] > 0) {
        // NumericVector temp = colcmm.row(ij);
        // Rcout << "The value of colcmm: " << temp << "\n";

        // temp = colc_cov_coef.row(cc);
        // Rcout << "The value of colc_cov_coef: " << temp << "\n";

        for (ll=0; ll < ncolccov; ll++) {
            linear_part += colcmm(ij,ll)*colc_cov_coef(cc,ll);
        }
        // Rcout << "The value of linear_part with colc_cov: " << linear_part << "\n";
    }
    if (paramlengths["cov"] > 0) {
        for (ll=0; ll < paramlengths["cov"]; ll++) {
            linear_part += covmm(ij,ll)*cov_coef[ll];
        }
        // Rcout << "The value of linear_part with cov: " << linear_part << "\n";
    }

    return linear_part;
}

double rcpp_theta_from_linear(const String model,
                              const double linear_part,
                              const int ymatij_idx,
                              const NumericVector mu,
                              const NumericVector phi,
                              const int & q,
                              const double & epsilon) {

    double theta = 0;
    NumericVector theta_all (q);
    double theta_sum = 0;

    int kk;

    theta_sum = 0;
    if (model == "OSM") {
        theta_all[0] = 1;
        theta_sum += 1;
        for (kk=1; kk < q; kk++) {
            theta_all[kk] = exp(mu[kk] + phi[kk]*linear_part);
            theta_sum += theta_all[kk];
        }
        theta = theta_all[ymatij_idx]/theta_sum;
    } else if (model == "POM") {
        theta_all[0] = rcpp_expit(mu[0] - linear_part);
        theta_sum = theta_all[0];
        for (kk=1; kk < q-1; kk++) {
            theta_all[kk] = rcpp_expit(mu[kk] - linear_part) -
                rcpp_expit(mu[kk-1] - linear_part);
            theta_sum += theta_all[kk];
        }
        theta_all[q-1] = 1-theta_sum;
        theta = theta_all[ymatij_idx];
    } else if (model == "Binary") {
        theta_all[0] = 1;
        theta_all[1] = exp(mu[0] + linear_part);
        theta_sum = 1 + theta_all[1];
        theta = theta_all[ymatij_idx]/theta_sum;
    }

    if (theta <= epsilon) {
        theta = epsilon;
    }

    return theta;
}

// [[Rcpp::export]]
double rcpp_Rclusterll(const NumericVector & invect,
                       const String & model,
                       const NumericMatrix & ydf,
                       const NumericMatrix & rowcmm,
                       const NumericMatrix & colcmm,
                       const NumericMatrix & covmm,
                       const NumericMatrix & pprm,
                       const NumericVector & piv,
                       const IntegerVector & paramlengths,
                       const int & RG, const int & p, const int & n, const int & q,
                       const double & epsilon, const bool & constraint_sum_zero,
                       const bool & partial, const bool & incomplete) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: paramlengths["rowc"], else you will get a
    // fatal error that will crash rr!

    // Note that invect MUST be the first argument, so this likelihood function
    // can be used inside optim()!

    NumericVector mu(q,NA_REAL);
    NumericVector phi(q,NA_REAL);

    NumericVector rowc_coef(RG);
    // Need to set up the objects even if we won't use them, because it's not
    // possible to set them up within if statements
    NumericVector col_coef(p,NA_REAL);

    NumericMatrix rowc_col_coef(RG,p);
    std::fill( rowc_col_coef.begin(), rowc_col_coef.end(), NA_REAL ) ;

    int nrowccov = paramlengths["rowc_cov"]/RG;
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;
    int ncolccov = 0;

    NumericVector cov_coef(paramlengths["cov"],NA_REAL);

    NumericVector colc_coef(1,NA_REAL);
    NumericVector row_coef(1,NA_REAL);
    NumericMatrix colc_row_coef(1,1);
    std::fill( colc_row_coef.begin(), colc_row_coef.end(), NA_REAL ) ;
    NumericMatrix colc_cov_coef(1,1);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;
    NumericMatrix rowc_colc_coef(1,1);
    std::fill( rowc_colc_coef.begin(), rowc_colc_coef.end(), NA_REAL ) ;
    int CG = 0;

    rcpp_unpack(model, invect, paramlengths, mu, phi, rowc_coef, colc_coef,
                rowc_colc_coef, row_coef, col_coef, rowc_col_coef, colc_row_coef,
                rowc_cov_coef, colc_cov_coef, cov_coef,
                RG, CG, p, n, q, constraint_sum_zero);

    // Rcout << "The value of mu : " << mu << "\n";
    // Rcout << "The value of phi: " << phi << "\n";
    // Rcout << "The value of rowc_coef : " << rowc_coef << "\n";
    // Rcout << "The value of col_coef : " << col_coef << "\n";
    // Rcout << "The value of rowc_col_coef : " << rowc_col_coef << "\n";
    // Rcout << "The value of rowc_cov_coef : " << rowc_cov_coef << "\n";
    // Rcout << "The value of cov_coef : " << cov_coef << "\n";
    //
    // Rcout << "The value of colc_coef : " << colc_coef << "\n";
    // Rcout << "The value of row_coef : " << row_coef << "\n";
    // Rcout << "The value of colc_row_coef : " << colc_row_coef << "\n";
    // Rcout << "The value of colc_cov_coef : " << colc_cov_coef << "\n";
    // Rcout << "The value of rowc_colc_coef : " << rowc_colc_coef << "\n";

    int rr, ij, ii, jj;
    int cc = 0;

    NumericVector yval;
    int ymatij_idx = 0;
    double linear_part = 0;
    double theta = 0;
    double log_thetaymat = 0;

    double logl;

    // Adjust very small values of piv to avoid errors when calculating log(piv)
    for (rr=0; rr < RG; rr++) {
        if (piv[rr] < epsilon) piv[rr] = epsilon;
    }

    if (!incomplete) {
        double llc = 0;

        for (rr=0; rr < RG; ++rr) {
            for (ij=0; ij < ydf.nrow(); ++ij) {

                ii = ydf(ij,1)-1;
                jj = ydf(ij,2)-1;

                // Rcout << "The R-based value of ii : " << ii+1 << "\n";
                // Rcout << "The R-based value of jj : " << jj+1 << "\n";
                yval = ydf(ij,0);
                if (all(is_finite(yval)) & all(!is_nan(yval))) {
                    linear_part = rcpp_linear_part(ydf, rowcmm, colcmm, covmm,
                                                   paramlengths,
                                                   rowc_coef, colc_coef,
                                                   rowc_colc_coef,
                                                   row_coef, col_coef,
                                                   rowc_col_coef, colc_row_coef,
                                                   rowc_cov_coef, colc_cov_coef,
                                                   cov_coef,
                                                   RG, CG, p, n, q,
                                                   nrowccov, ncolccov,
                                                   rr, cc, ij, ii, jj);

                    ymatij_idx = ydf(ij,0)-1;

                    theta = rcpp_theta_from_linear(model, linear_part, ymatij_idx, mu, phi, q, epsilon);

                    log_thetaymat = log(theta);
                } else {
                    log_thetaymat = 0;
                }
                Rcout << "The value of log_thetaymat : " << log_thetaymat << "\n";

                // Rcout << "The value of llc component : " << pprm(ii,rr)*log_thetaymat << "\n";
                llc += pprm(ii,rr)*log_thetaymat;
            }
        }

        if (!partial) {
            llc += sum(pprm*log(piv));
        }

        logl = llc;
    } else {
        logl = 0;

        NumericVector log_components (RG);
        double log_sumoverR = 0;
        double theta;
        double log_theta;

        for (ii=0; ii < n; ++ii) {
            Rcout << "The R-based value of ii : " << ii+1 << "\n";

            log_components.fill(0);

            for (rr=0; rr < RG; rr++) {
                log_components[rr] = log(piv[rr]);

                for (ij=0; ij < ydf.nrow(); ij++) {
                    if (ydf(ij,1) == ii+1) {

                        jj = ydf(ij,2)-1;
                        Rcout << "The R-based value of jj : " << jj+1 << "\n";
                        yval = ydf(ij,0);
                        if (all(is_finite(yval)) & all(!is_nan(yval))) {
                            linear_part = rcpp_linear_part(ydf, rowcmm, colcmm, covmm,
                                                           paramlengths,
                                                           rowc_coef, colc_coef,
                                                           rowc_colc_coef,
                                                           row_coef, col_coef,
                                                           rowc_col_coef, colc_row_coef,
                                                           rowc_cov_coef, colc_cov_coef,
                                                           cov_coef,
                                                           RG, CG, p, n, q,
                                                           nrowccov, ncolccov,
                                                           rr, cc, ij, ii, jj);
                            ymatij_idx = ydf(ij,0)-1;

                            theta = rcpp_theta_from_linear(model, linear_part, ymatij_idx, mu, phi, q);
                            log_theta = log(theta);
                            if (!NumericVector::is_na(log_theta) &&
                                !Rcpp::traits::is_nan<REALSXP>(log_theta) &&
                                !Rcpp::traits::is_infinite<REALSXP>(log_theta)) {
                                log_components[rr] += log_theta;
                            }
                        }
                    }
                }
            }
            log_sumoverR = log(sum(exp(log_components - max(log_components)))) + max(log_components);
            logl += log_sumoverR;
        }

        if (logl == 0) logl = -1E-40;
    }

    return logl;
}

// [[Rcpp::export]]
double rcpp_Biclusterll(const NumericVector & invect,
                        const String & model,
                        const NumericMatrix & ydf,
                        const NumericMatrix & rowcmm,
                        const NumericMatrix & colcmm,
                        const NumericMatrix & covmm,
                        const NumericMatrix & pprm,
                        const NumericMatrix & ppcm,
                        const NumericVector & piv,
                        const NumericVector & kappav,
                        const IntegerVector & paramlengths,
                        const int & RG, const int & CG,
                        const int & p, const int & n, const int & q,
                        const double & epsilon,
                        const bool & constraint_sum_zero, const bool & partial,
                        const bool & incomplete, double & llc) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: paramlengths["rowc"], else you will get a
    // fatal error that will crash rr!

    // Note that invect MUST be the first argument, so this likelihood function
    // can be used inside optim()!

    NumericVector mu(q,NA_REAL);
    NumericVector phi(q,NA_REAL);

    NumericVector rowc_coef(RG);
    NumericVector colc_coef(CG);

    NumericMatrix rowc_colc_coef(RG,CG);
    std::fill( rowc_colc_coef.begin(), rowc_colc_coef.end(), NA_REAL ) ;

    // Need to set up the objects even if we won't use them, because it's not
    // possible to set them up within if statements
    NumericVector col_coef(p,NA_REAL);
    NumericVector row_coef(n,NA_REAL);

    NumericMatrix rowc_col_coef(RG,p);
    std::fill( rowc_col_coef.begin(), rowc_col_coef.end(), NA_REAL ) ;

    NumericMatrix colc_row_coef(CG,n);
    std::fill( colc_row_coef.begin(), colc_row_coef.end(), NA_REAL ) ;

    int nrowccov = 0;
    if (paramlengths["rowc_cov"] > 0) {
        nrowccov = paramlengths["rowc_cov"]/RG;
    }
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;

    int ncolccov = 0;
    if (paramlengths["colc_cov"] > 0) {
        ncolccov = paramlengths["colc_cov"]/CG;
    }
    // Rcout << "The value of ncolccov: " << ncolccov << "\n";
    NumericMatrix colc_cov_coef(CG,ncolccov);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;

    NumericVector cov_coef(paramlengths["cov"],NA_REAL);

    rcpp_unpack(model, invect, paramlengths, mu, phi, rowc_coef, colc_coef,
                rowc_colc_coef, row_coef, col_coef, rowc_col_coef, colc_row_coef,
                rowc_cov_coef, colc_cov_coef, cov_coef,
                RG, CG, p, n, q, constraint_sum_zero);

    // Rcout << "The value of mu : " << mu << "\n";
    // Rcout << "The value of phi: " << phi << "\n";
    // Rcout << "The value of rowc_coef : " << rowc_coef << "\n";
    // Rcout << "The value of col_coef : " << col_coef << "\n";
    // Rcout << "The value of rowc_col_coef : " << rowc_col_coef << "\n";
    // Rcout << "The value of rowc_cov_coef : " << rowc_cov_coef << "\n";
    // Rcout << "The value of cov_coef : " << cov_coef << "\n";

    // Rcout << "The value of colc_coef : " << colc_coef << "\n";
    // Rcout << "The value of row_coef : " << row_coef << "\n";
    // Rcout << "The value of colc_row_coef : " << colc_row_coef << "\n";
    // Rcout << "The value of colc_cov_coef : " << colc_cov_coef << "\n";
    // Rcout << "The value of rowc_colc_coef : " << rowc_colc_coef << "\n";

    int rr, cc, ij, ii, jj;

    NumericVector yval;
    int ymatij_idx = 0;
    double linear_part = 0;
    double theta = 0;
    double log_thetaymat;

    double logl = 0;

    // Adjust very small values of piv to avoid errors when calculating log(piv)
    for (rr=0; rr < RG; rr++) {
        if (piv[rr] < epsilon) piv[rr] = epsilon;
    }
    // Similarly for kappav
    for (cc=0; cc < CG; cc++) {
        if (kappav[cc] < epsilon) kappav[cc] = epsilon;
    }

    if (!incomplete || NumericVector::is_na(llc)) {
        llc = 0;

        for (rr=0; rr < RG; ++rr) {
            for (cc=0; cc < CG; ++cc) {

                for (ij=0; ij < ydf.nrow(); ++ij) {
                    ii = ydf(ij,1)-1;
                    jj = ydf(ij,2)-1;
                    // Rcout << "The R-based value of ii : " << ii+1 << "\n";
                    // Rcout << "The R-based value of jj : " << jj+1 << "\n";

                    yval = ydf(ij,0);
                    if (all(is_finite(yval)) & all(!is_nan(yval))) {
                        ymatij_idx = ydf(ij,0)-1;

                        // Rcout << "The value of ymatij_idx : " << ymatij_idx << "\n";

                        linear_part = rcpp_linear_part(ydf, rowcmm, colcmm, covmm,
                                                       paramlengths,
                                                       rowc_coef, colc_coef,
                                                       rowc_colc_coef,
                                                       row_coef, col_coef,
                                                       rowc_col_coef, colc_row_coef,
                                                       rowc_cov_coef, colc_cov_coef,
                                                       cov_coef,
                                                       RG, CG, p, n, q,
                                                       nrowccov, ncolccov,
                                                       rr, cc, ij, ii, jj);

                        theta = rcpp_theta_from_linear(model, linear_part, ymatij_idx, mu, phi, q);
                        // Rcout << "The value of theta : " << theta << "\n";

                        log_thetaymat = log(theta);
                    } else {
                        log_thetaymat = 0;
                    }
                    // Rcout << "The value of log_thetaymat : " << log_thetaymat << "\n";
                    // Rcout << "The value of llc component : " << pprm(ii,rr)*log_thetaymat*ppcm(jj,cc) << "\n";
                    llc += pprm(ii,rr)*log_thetaymat*ppcm(jj,cc);
                }
            }
        }

        if (!partial) {
            llc += sum(pprm*log(piv));
            llc += sum(ppcm*log(kappav));
        }

        logl = llc;
    }

    if (incomplete) {
        double logl_correction = 0;

        NumericMatrix tau_numerator (RG,CG);
        double this_num;
        double tau_denominator;
        double logl_correction_term;

        for (ij=0; ij < ydf.nrow(); ++ij) {
            ii = ydf(ij,1)-1;
            jj = ydf(ij,2)-1;

            yval = ydf(ij,0);
            if (all(is_finite(yval)) & all(!is_nan(yval))) {
                ymatij_idx = ydf(ij,0)-1;

                tau_denominator = 0;

                for (rr=0; rr < RG; rr++) {
                    for (cc=0; cc < CG; cc++) {
                        linear_part = rcpp_linear_part(ydf, rowcmm, colcmm, covmm,
                                                       paramlengths,
                                                       rowc_coef, colc_coef,
                                                       rowc_colc_coef,
                                                       row_coef, col_coef,
                                                       rowc_col_coef, colc_row_coef,
                                                       rowc_cov_coef, colc_cov_coef,
                                                       cov_coef,
                                                       RG, CG, p, n, q,
                                                       nrowccov, ncolccov,
                                                       rr, cc, ij, ii, jj);

                        theta = rcpp_theta_from_linear(model, linear_part, ymatij_idx, mu, phi, q, epsilon);

                        tau_numerator(rr,cc) = piv[rr]*theta*kappav[cc];
                        tau_denominator += tau_numerator(rr,cc);
                    }
                }

                for (rr=0; rr < RG; rr++) {
                    for (cc=0; cc < CG; cc++) {

                        this_num = pprm(ii,rr);
                        Rcout << "This z_ir is:" << this_num << "\n";
                        this_num = ppcm(jj,cc);
                        Rcout << "This x_jc is:" << this_num << "\n";

                        logl_correction_term = pprm(ii,rr)*ppcm(jj,cc)*log(tau_numerator(rr,cc)/tau_denominator);
                        if (!NumericVector::is_na(logl_correction_term) &&
                            !Rcpp::traits::is_nan<REALSXP>(logl_correction_term) &&
                            !Rcpp::traits::is_infinite<REALSXP>(logl_correction_term)) {
                            logl_correction += logl_correction_term;
                            Rcout << "The logl correction term is:" << logl_correction_term << "\n";
                        }
                    }
                }
            }
        }

        if (!NumericVector::is_na(logl_correction) &&
            !Rcpp::traits::is_nan<REALSXP>(logl_correction) &&
            !Rcpp::traits::is_infinite<REALSXP>(logl_correction)) {
            logl -= logl_correction;
        } else {
            logl = R_NegInf;
        }
    }

    return logl;
}

// [[Rcpp::export]]
NumericMatrix rcpp_Rcluster_Estep(const NumericVector & invect,
                                  const String & model,
                                  const NumericMatrix & ydf,
                                  const NumericMatrix & rowcmm,
                                  const NumericMatrix & colcmm,
                                  const NumericMatrix & covmm,
                                  const NumericVector & piv,
                                  const IntegerVector & paramlengths,
                                  const int & RG, const int & p, const int & n, const int & q,
                                  const double & epsilon,
                                  const bool & constraint_sum_zero) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: paramlengths["rowc"], else you will get a
    // fatal error that will crash rr!

    // Note that invect MUST be the first argument, so this likelihood function
    // can be used inside optim()!

    NumericVector mu(q,NA_REAL);
    NumericVector phi(q,NA_REAL);

    NumericVector rowc_coef(RG);
    // Need to set up the objects even if we won't use them, because it's not
    // possible to set them up within if statements
    NumericVector col_coef(p,NA_REAL);

    NumericMatrix rowc_col_coef(RG,p);
    std::fill( rowc_col_coef.begin(), rowc_col_coef.end(), NA_REAL ) ;

    int nrowccov = paramlengths["rowc_cov"]/RG;
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;
    int ncolccov = 0;

    NumericVector cov_coef(paramlengths["cov"],NA_REAL);

    NumericVector colc_coef(1,NA_REAL);
    NumericVector row_coef(1,NA_REAL);
    NumericMatrix colc_row_coef(1,1);
    std::fill( colc_row_coef.begin(), colc_row_coef.end(), NA_REAL ) ;
    NumericMatrix colc_cov_coef(1,1);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;
    NumericMatrix rowc_colc_coef(1,1);
    std::fill( rowc_colc_coef.begin(), rowc_colc_coef.end(), NA_REAL ) ;
    int CG = 0;

    rcpp_unpack(model, invect, paramlengths, mu, phi, rowc_coef, colc_coef,
                rowc_colc_coef, row_coef, col_coef, rowc_col_coef, colc_row_coef,
                rowc_cov_coef, colc_cov_coef, cov_coef,
                RG, CG, p, n, q, constraint_sum_zero);

    int rr, ij, ii, jj;
    int cc=0;

    NumericMatrix ppm(n, RG);
    NumericMatrix ppm_raw(n, RG);
    NumericVector ppm_row(RG);
    NumericVector ppm_row_adjusted(RG);
    double ppm_rowminabs;
    double ppm_log_sum_exp_adjusted;

    IntegerMatrix ppm_started(n,RG);

    NumericVector yval;
    int ymatij_idx;
    double linear_part;
    double theta;
    double log_thetaymat;

    // Adjust very small values of piv to avoid errors when calculating log(piv)
    for (rr=0; rr < RG; rr++) {
        if (piv[rr] < epsilon) piv[rr] = epsilon;
    }

    for (ij=0; ij < ydf.nrow(); ij++) {

        ii = ydf(ij,1)-1;
        jj = ydf(ij,2)-1;

        for (rr=0; rr < RG; rr++) {
            if (ppm_started(ii,rr) == 0) {
                ppm_raw(ii,rr) = log(piv[rr]);
                ppm_started(ii,rr) = 1;
            }

            yval = ydf(ij,0);
            if (all(is_finite(yval)) & all(!is_nan(yval))) {
                linear_part = rcpp_linear_part(ydf, rowcmm, colcmm, covmm,
                                               paramlengths,
                                               rowc_coef, colc_coef,
                                               rowc_colc_coef,
                                               row_coef, col_coef,
                                               rowc_col_coef, colc_row_coef,
                                               rowc_cov_coef, colc_cov_coef,
                                               cov_coef,
                                               RG, CG, p, n, q,
                                               nrowccov, ncolccov,
                                               rr, cc, ij, ii, jj);

                ymatij_idx = ydf(ij,0)-1;

                theta = rcpp_theta_from_linear(model, linear_part, ymatij_idx, mu, phi, q, epsilon);

                log_thetaymat = log(theta);
            } else {
                log_thetaymat = 0;
            }
            ppm_raw(ii,rr) += log_thetaymat;
        }
    }
    for (ii=0; ii < n; ii++) {
        ppm_row = ppm_raw.row(ii);
        ppm_rowminabs = min(abs(ppm_row));
        ppm_row_adjusted = ppm_row + ppm_rowminabs;
        ppm_log_sum_exp_adjusted = log(sum(exp(ppm_row_adjusted)));
        ppm_row = ppm_row - ppm_log_sum_exp_adjusted + ppm_rowminabs;
        ppm_row = exp(ppm_row);
        ppm.row(ii) = ppm_row/sum(ppm_row);
    }
    // Rcout << "The value of ppm : " << ppm << "\n";

    return ppm;
}

// [[Rcpp::export]]
NumericMatrix rcpp_Bicluster_Estep(const NumericVector & invect,
                                   const String & model,
                                   const NumericMatrix & ydf,
                                   const NumericMatrix & rowcmm,
                                   const NumericMatrix & colcmm,
                                   const NumericMatrix & covmm,
                                   const NumericVector & piv,
                                   const NumericVector & kappav,
                                   const IntegerVector & paramlengths,
                                   const int & RG, const int & CG,
                                   const int & p, const int & n, const int & q,
                                   const double & epsilon,
                                   const bool & constraint_sum_zero,
                                   const bool & row_clusters) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: paramlengths["rowc"], else you will get a
    // fatal error that will crash rr!

    // Note that invect MUST be the first argument, so this likelihood function
    // can be used inside optim()!

    NumericVector mu(q,NA_REAL);
    NumericVector phi(q,NA_REAL);

    NumericVector rowc_coef(RG);
    NumericVector colc_coef(CG);

    NumericMatrix rowc_colc_coef(RG,CG);
    std::fill( rowc_colc_coef.begin(), rowc_colc_coef.end(), NA_REAL ) ;

    // Need to set up the objects even if we won't use them, because it's not
    // possible to set them up within if statements
    NumericVector col_coef(p,NA_REAL);
    NumericVector row_coef(n,NA_REAL);

    NumericMatrix rowc_col_coef(RG,p);
    std::fill( rowc_col_coef.begin(), rowc_col_coef.end(), NA_REAL ) ;

    NumericMatrix colc_row_coef(CG,n);
    std::fill( colc_row_coef.begin(), colc_row_coef.end(), NA_REAL ) ;

    int nrowccov = 0;
    if (paramlengths["rowc_cov"] > 0) {
        nrowccov = paramlengths["rowc_cov"]/RG;
    }
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;

    int ncolccov = 0;
    if (paramlengths["colc_cov"] > 0) {
        ncolccov = paramlengths["colc_cov"]/CG;
    }
    // Rcout << "The value of ncolccov: " << ncolccov << "\n";
    NumericMatrix colc_cov_coef(CG,ncolccov);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;

    NumericVector cov_coef(paramlengths["cov"],NA_REAL);

    rcpp_unpack(model, invect, paramlengths, mu, phi, rowc_coef, colc_coef,
                rowc_colc_coef, row_coef, col_coef, rowc_col_coef, colc_row_coef,
                rowc_cov_coef, colc_cov_coef, cov_coef,
                RG, CG, p, n, q, constraint_sum_zero);

    int rr, cc, ij, ii, jj, ll;

    NumericVector yval;
    int ymatij_idx = 0;
    double linear_part = 0;
    double theta = 0;
    double sum_pikappa_theta;

    int nelements;
    int nclust;
    if (row_clusters) {
        nelements = n;
        nclust = RG;
    } else {
        nelements = p;
        nclust = CG;
    }
    // Rcout << "The value of nelements : " << nelements << "\n";
    // Rcout << "The value of nclust : " << nclust << "\n";

    NumericMatrix ppm(nelements, nclust);
    NumericMatrix ppm_raw(nelements, nclust);
    NumericVector ppm_row(nclust);
    NumericVector ppm_row_adjusted(nclust);
    double ppm_rowminabs;
    double ppm_log_sum_exp_adjusted;

    NumericMatrix ppm_started(nelements, nclust);

    // Adjust very small values of piv to avoid errors when calculating log(piv)
    for (rr=0; rr < RG; rr++) {
        if (piv[rr] < epsilon) piv[rr] = epsilon;
    }
    // Similarly for kappav
    for (cc=0; cc < CG; cc++) {
        if (kappav[cc] < epsilon) kappav[cc] = epsilon;
    }

    for (ij=0; ij < ydf.nrow(); ij++) {

        ii = ydf(ij,1)-1;
        jj = ydf(ij,2)-1;
        // Rcout << "The R-based value of ii : " << ii+1 << "\n";
        // Rcout << "The R-based value of jj : " << jj+1 << "\n";

        if (row_clusters) {

            for (rr=0; rr < RG; rr++) {
                // Need this step because we are looping over ij, but for each
                // ii and rr we only want to include one log(piv[rr]), instead of
                // adding one for each jj
                if (ppm_started(ii,rr) == 0) {
                    ppm_raw(ii,rr) = log(piv[rr]);
                    ppm_started(ii,rr) = 1;
                }

                yval = ydf(ij,0);
                if (all(is_finite(yval)) & all(!is_nan(yval))) {
                    ymatij_idx = ydf(ij,0)-1;
                    sum_pikappa_theta = 0;

                    for (cc=0; cc < CG; cc++) {
                        linear_part = rcpp_linear_part(ydf, rowcmm, colcmm, covmm,
                                                       paramlengths,
                                                       rowc_coef, colc_coef,
                                                       rowc_colc_coef,
                                                       row_coef, col_coef,
                                                       rowc_col_coef, colc_row_coef,
                                                       rowc_cov_coef, colc_cov_coef,
                                                       cov_coef,
                                                       RG, CG, p, n, q,
                                                       nrowccov, ncolccov,
                                                       rr, cc, ij, ii, jj);

                        theta = rcpp_theta_from_linear(model, linear_part, ymatij_idx, mu, phi, q, epsilon);

                        sum_pikappa_theta += kappav[cc]*theta;
                    }

                    // Rcout << "The value of the additional log term : " << log(sum_pikappa_theta) << "\n";
                    ppm_raw(ii,rr) += log(sum_pikappa_theta);
                }
            }
        } else {
            for (cc=0; cc < CG; cc++) {
                // Need this step because we are looping over ij, but for each
                // jj and cc we only want to include one log(kappav[cc]),
                // instead of adding one for each ii
                if (ppm_started(jj,cc) == 0) {
                    ppm_raw(jj,cc) = log(kappav[cc]);
                    ppm_started(jj,cc) = 1;
                }

                yval = ydf(ij,0);
                if (all(is_finite(yval)) & all(!is_nan(yval))) {
                    ymatij_idx = ydf(ij,0)-1;
                    sum_pikappa_theta = 0;

                    for (rr=0; rr < RG; rr++) {
                        linear_part = rcpp_linear_part(ydf, rowcmm, colcmm, covmm,
                                                       paramlengths,
                                                       rowc_coef, colc_coef,
                                                       rowc_colc_coef,
                                                       row_coef, col_coef,
                                                       rowc_col_coef, colc_row_coef,
                                                       rowc_cov_coef, colc_cov_coef,
                                                       cov_coef,
                                                       RG, CG, p, n, q,
                                                       nrowccov, ncolccov,
                                                       rr, cc, ij, ii, jj);

                        theta = rcpp_theta_from_linear(model, linear_part, ymatij_idx, mu, phi, q, epsilon);

                        sum_pikappa_theta += piv[rr]*theta;
                    }
                    ppm_raw(jj,cc) += log(sum_pikappa_theta);
                }
            }
        }
    }
    // Rcout << "The value of ppm_raw : " << ppm_raw << "\n";

    for (ll=0; ll < nelements; ll++) {
        ppm_row = ppm_raw.row(ll);
        ppm_rowminabs = min(abs(ppm_row));
        ppm_row_adjusted = ppm_row + ppm_rowminabs;
        ppm_log_sum_exp_adjusted = log(sum(exp(ppm_row_adjusted)));
        ppm_row = ppm_row - ppm_log_sum_exp_adjusted + ppm_rowminabs;
        ppm_row = exp(ppm_row);
        ppm.row(ll) = ppm_row/sum(ppm_row);
    }

    return ppm;
}


/*** R

*/


