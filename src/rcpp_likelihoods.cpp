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

void rcpp_unpack(const int & model_num,
                 const NumericVector & invect,
                 const IntegerVector & param_lengths,
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
    // vector, only double quotes: param_lengths["rowc"], else you will get a
    // fatal error that will crash R;

    // param_lengths entries: c('mu','phi','rowc','colc','rowc_colc','row','col','rowc_col','colc_row','rowc_cov','colc_cov','cov')]

    int ind, kk, rr, cc, ii, jj, ll;
    int nelts = 0;

    if (model_num == 1) {
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

    } else if (model_num == 2) {
        // Convert to mu from w, where w can vary between -Inf and +Inf
        // but mu must be increasing i.e. mu[1] <= mu[2] <= mu[3]...
        mu[0] = invect[0];
        for (kk=1; kk < q-1; kk++) {
            mu[kk] = mu[kk-1] + exp(invect[kk]);
        }
        nelts += q-1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of mu : " << mu << "\n";
    } else if (model_num == 3) {
        mu[0] = invect[0];
        nelts += 1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of mu : " << mu << "\n";
    }

    if (param_lengths[2] > 0) {
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
    if (param_lengths[3] > 0) {
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
    if (param_lengths[4] > 0) {

        if (param_lengths[2] > 0 && param_lengths[3] > 0) {
            std::fill( rowc_colc_coef.row(RG-1).begin(), rowc_colc_coef.row(RG-1).end(), 0 );
            std::fill( rowc_colc_coef.column(CG-1).begin(), rowc_colc_coef.column(CG-1).end(), 0 );

            for (rr=0; rr < RG-1; rr++) {
                for (cc=0; cc < CG-1; cc++) {
                    ind = nelts + rr*(CG-1)+cc;

                    // Rcout << "The index of invect : " << ind << "\n";
                    // Using constraint formulation from original POM code, with
                    // final row of rowc.colc.coef equal to negative sum of other
                    // rows. This is unlike the v0.1 clustord code and the original
                    // OSM code, had FIRST row of rowc.colc.coef equal to negative
                    // sum of other rows
                    rowc_colc_coef(rr,cc) = invect[ind];
                    rowc_colc_coef(rr,CG-1) -= invect[ind]; // fill last column with negative sum of other columns
                    rowc_colc_coef(RG-1,cc) -= invect[ind]; // fill last row with negative sum of other rows
                    rowc_colc_coef(RG-1,CG-1) += invect[ind]; // fill last element with sum of all independent elements
                }
            }
            nelts += (RG-1)*(CG-1);
        } else {
            rowc_colc_coef(RG-1,CG-1) = 0;

            for (rr=0; rr < RG; rr++) {
                for (cc=0; cc < CG; cc++) {

                    if (rr != RG-1 || cc != CG-1) {
                        ind = nelts + rr*CG + cc;
                        // Rcout << "The index of invect : " << ind << "\n";
                        rowc_colc_coef(rr,cc) = invect[ind];
                        // fill last element with negative sum of other elements
                        rowc_colc_coef(RG-1,CG-1) -= invect[ind];
                    }
                }
            }
            nelts += (RG*CG-1);
        }
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of rowc_colc_coef : " << rowc_colc_coef << "\n";
    }

    if (param_lengths[5] > 0) {
        if (constraint_sum_zero) {
            row_coef[n-1] = 0;
            for (ii=0; ii < n-1; ii++) {
                ind = nelts + ii;
                // Rcout << "The index of invect : " << ind << "\n";
                row_coef[ii] = invect[ind];
                row_coef[n-1] -= invect[ind];
            }
        } else {
            row_coef[0] = 0;
            for (ii=0; ii < n-1; ii++) {
                ind = nelts + ii;
                row_coef[ii+1] = invect[ind];
            }
        }
        nelts += n-1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of row_coef : " << row_coef << "\n";
    }
    if (param_lengths[6] > 0) {
        if (constraint_sum_zero) {
            col_coef[p-1] = 0;
            for (jj=0; jj < p-1; jj++) {
                ind = nelts + jj;
                // Rcout << "The index of invect : " << ind << "\n";
                col_coef[jj] = invect[ind];
                col_coef[p-1] -= invect[ind];
                // Rcout << "The value of col_coef : " << col_coef << "\n";
            }
        } else {
            col_coef[0] = 0;
            for (jj=0; jj < p-1; jj++) {
                ind = nelts + jj;
                col_coef[jj+1] = invect[ind];
            }
        }
        nelts += p-1;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of col_coef : " << col_coef << "\n";
    }

    if (param_lengths[7] > 0) {
        // NOTE: have to make sure to fill the matrix of rowc_col coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        if (param_lengths[2] > 0) {
            std::fill( rowc_col_coef.row(RG-1).begin(), rowc_col_coef.row(RG-1).end(), 0 );
            std::fill( rowc_col_coef.column(p-1).begin(), rowc_col_coef.column(p-1).end(), 0 );

            for (rr=0; rr < RG-1; rr++) {
                for (jj=0; jj < p-1; jj++) {
                    ind = nelts + rr*(p-1) + jj;
                    // Rcout << "The index of invect : " << ind << "\n";

                    rowc_col_coef(rr,jj) = invect[ind];
                    rowc_col_coef(rr,p-1) -= invect[ind]; // fill last column with negative sum of other columns
                    rowc_col_coef(RG-1,jj) -= invect[ind]; // fill last row with negative sum of other rows
                    rowc_col_coef(RG-1,p-1) += invect[ind]; // fill last element with sum of all independent elements
                }
            }
            nelts += (RG-1)*(p-1);
        } else {
            rowc_col_coef(RG-1,p-1) = 0;
            for (rr=0; rr < RG; rr++) {
                for (jj=0; jj < p; jj++) {
                    if (rr != RG-1 || jj != p-1) {
                        ind = nelts + rr*p + jj;
                        // Rcout << "The index of invect : " << ind << "\n";
                        rowc_col_coef(rr,jj) = invect[ind];
                        rowc_col_coef(RG-1,p-1) -= invect[ind];
                    }
                }
            }
            nelts += RG*p-1;
        }
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of rowc_col_coef : " << rowc_col_coef << "\n";
    }
    if (param_lengths[8] > 0) {
        // NOTE: have to make sure to fill the matrix of colc_row coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        colc_row_coef(CG-1,n-1) = 0;
        if (param_lengths[3] > 0) {
            std::fill( rowc_colc_coef.row(CG-1).begin(), rowc_colc_coef.row(CG-1).end(), 0 );
            std::fill( rowc_colc_coef.column(n-1).begin(), rowc_colc_coef.column(n-1).end(), 0 );

            for (cc=0; cc < CG-1; cc++) {
                for (ii=0; ii < n-1; ii++) {
                    ind = nelts + cc*(n-1) + ii;

                    colc_row_coef(cc,ii) = invect[ind];
                    colc_row_coef(cc,n-1) -= invect[ind]; // fill last column with negative sum of other columns
                    colc_row_coef(CG-1,ii) -= invect[ind]; // fill last row with negative sum of other rows
                    colc_row_coef(CG-1,n-1) += invect[ind]; // fill last element with sum of all independent elements
                }
            }
            nelts += (CG-1)*(n-1);
        } else {
            colc_row_coef(CG-1,n-1) = 0;
            for (cc=0; cc < CG; cc++) {
                for (ii=0; ii < n; ii++) {
                    if (cc != CG-1 || ii != n-1) {
                        ind = nelts + cc*n + ii;
                        colc_row_coef(cc,ii) = invect[ind];
                        colc_row_coef(CG-1,n-1) -= invect[ind];
                    }
                }
            }
            nelts += CG*n-1;
        }
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of colc_row_coef : " << colc_row_coef << "\n";
    }

    if (param_lengths[9] > 0 && RG > 0) {
        // NOTE: have to make sure to fill the matrix of rowc_cov coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        int nrowccov = param_lengths[9]/RG;
        // Rcout << "The value of nrowccov : " << nrowccov << "\n";
        for (rr=0; rr < RG; rr++) {
            for (ll=0; ll < nrowccov; ll++) {
                ind = nelts + rr*nrowccov + ll;
                // Rcout << "The index of invect : " << ind << "\n";
                rowc_cov_coef(rr,ll) = invect[ind];
            }
        }
        nelts += RG*nrowccov;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of rowc_cov_coef : " << rowc_cov_coef << "\n";
    }
    if (param_lengths[10] > 0 && CG > 0) {
        // NOTE: have to make sure to fill the matrix of colc_cov coefs from the
        // parameter vector IN THE SAME WAY as the matrix is filled in the rr
        // unpack_parvec function!
        int ncolccov = param_lengths[10]/CG;
        for (cc=0; cc < CG; cc++) {
            for (ll=0; ll < ncolccov; ll++) {
                ind = nelts + cc*ncolccov + ll;
                // Rcout << "The index of invect : " << ind << "\n";
                colc_cov_coef(cc,ll) = invect[ind];
            }
        }
        nelts += CG*ncolccov;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of colc_cov_coef : " << colc_cov_coef << "\n";
    }

    int ncov = param_lengths[11];
    if (ncov > 0) {
        for (ll=0; ll < ncov; ll++) {
            ind = nelts + ll;
            // Rcout << "The index of invect : " << ind << "\n";
            cov_coef[ll] = invect[ind];
        }
        nelts += ncov;
        // Rcout << "The value of nelts : " << nelts << "\n";
        // Rcout << "The value of cov_coef : " << cov_coef << "\n";
    }
}

double rcpp_linear_part(const NumericMatrix ydf,
                        const NumericMatrix & rowc_mm,
                        const NumericMatrix & colc_mm,
                        const NumericMatrix & cov_mm,
                        const IntegerVector & param_lengths,
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

    double temp;

    // param_lengths entries: c('mu','phi','rowc','colc','rowc_colc','row','col','rowc_col','colc_row','rowc_cov','colc_cov','cov')]
    // Rcout << "The R-based value of ii : " << ii+1 << "\n";
    // Rcout << "The R-based value of jj : " << jj+1 << "\n";

    if (param_lengths[2] > 0) {
        linear_part += rowc_coef[rr];
        // Rcout << "The value of linear_part with rowc: " << linear_part << "\n";
    }
    if (param_lengths[3] > 0) {
        linear_part += colc_coef[cc];
        // Rcout << "The value of linear_part with colc: " << linear_part << "\n";
    }
    if (param_lengths[4] > 0) {
        linear_part += rowc_colc_coef(rr,cc);
        // Rcout << "The value of linear_part with rowc_colc: " << linear_part << "\n";
    }

    if (param_lengths[5] > 0) {
        linear_part += row_coef[ii];
        // Rcout << "The value of linear_part with row: " << linear_part << "\n";
    }
    if (param_lengths[6] > 0) {
        linear_part += col_coef[jj];
        // Rcout << "The value of linear_part with col: " << linear_part << "\n";
    }
    if (param_lengths[7] > 0) {
        linear_part += rowc_col_coef(rr,jj);
        // Rcout << "The value of linear_part with rowc_col: " << linear_part << "\n";
    }
    if (param_lengths[8] > 0) {
        linear_part += colc_row_coef(cc,ii);
        // Rcout << "The value of linear_part with colc_row: " << linear_part << "\n";
    }

    if (param_lengths[9] > 0) {
        // NumericVector temp = rowc_mm.row(ij);
        // Rcout << "The value of rowc_mm: " << temp << "\n";

        for (ll=0; ll < nrowccov; ll++) {
            linear_part += rowc_mm(ij,ll)*rowc_cov_coef(rr,ll);
        }
        // Rcout << "The value of linear_part with rowc_cov: " << linear_part << "\n";
    }
    if (param_lengths[10] > 0) {
        // NumericVector temp = colc_mm.row(ij);
        // Rcout << "The value of colc_mm: " << temp << "\n";

        // temp = colc_cov_coef.row(cc);
        // Rcout << "The value of colc_cov_coef: " << temp << "\n";

        for (ll=0; ll < ncolccov; ll++) {
            linear_part += colc_mm(ij,ll)*colc_cov_coef(cc,ll);
            if (std::isnan(linear_part)) {
                temp = colc_mm(ij,ll);
                Rcout << "The value of colc_mm element: " << temp << "\n";
                temp = colc_cov_coef(cc,ll);
                Rcout << "The value of colc_cov_coef element: " << temp << "\n";
            }
        }
        // Rcout << "The value of linear_part with colc_cov: " << linear_part << "\n";
    }
    if (param_lengths[11] > 0) {
        for (ll=0; ll < param_lengths[11]; ll++) {
            linear_part += cov_mm(ij,ll)*cov_coef[ll];
        }
        // Rcout << "The value of linear_part with cov: " << linear_part << "\n";
    }

    return linear_part;
}

double rcpp_theta_from_linear(const int & model_num,
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
    if (model_num == 1) {
        theta_all[0] = 1;
        theta_sum += 1;
        for (kk=1; kk < q; kk++) {
            theta_all[kk] = exp(mu[kk] + phi[kk]*linear_part);
            if (!std::isfinite(theta_all[kk])) {
                theta_all[kk] = 1;
            }

            theta_sum += theta_all[kk];
        }
        theta = theta_all[ymatij_idx]/theta_sum;
    } else if (model_num == 2) {
        theta_all[0] = rcpp_expit(mu[0] - linear_part);
        theta_sum = theta_all[0];
        for (kk=1; kk < q-1; kk++) {
            theta_all[kk] = rcpp_expit(mu[kk] - linear_part) -
                rcpp_expit(mu[kk-1] - linear_part);
            theta_sum += theta_all[kk];
        }
        theta_all[q-1] = 1-theta_sum;
        theta = theta_all[ymatij_idx];
    } else if (model_num == 3) {
        theta_all[0] = 1;
        theta_all[1] = exp(mu[0] + linear_part);
        theta_sum = 1 + theta_all[1];
        theta = theta_all[ymatij_idx]/theta_sum;
    }

    if (theta <= epsilon) {
        theta = epsilon;
    }

    if (std::isnan(theta)) {
        Rcout << "theta nan - The value of linear_part : " << linear_part << "\n";
    }

    return theta;
}

// [[Rcpp::export]]
double rcpp_Rclusterll(const NumericVector & invect,
                       const int & model_num,
                       const NumericMatrix & ydf,
                       const NumericMatrix & rowc_mm,
                       const NumericMatrix & colc_mm,
                       const NumericMatrix & cov_mm,
                       const NumericMatrix & ppr_m,
                       NumericVector pi_v,
                       const IntegerVector & param_lengths,
                       const int & RG, const int & p, const int & n, const int & q,
                       const double & epsilon, const bool & constraint_sum_zero,
                       const bool & partial, const bool & incomplete) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: param_lengths["rowc"], else you will get a
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

    int nrowccov = param_lengths["rowc_cov"]/RG;
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;
    int ncolccov = 0;

    NumericVector cov_coef(param_lengths["cov"],NA_REAL);

    NumericVector colc_coef(1,NA_REAL);
    NumericVector row_coef(1,NA_REAL);
    NumericMatrix colc_row_coef(1,1);
    std::fill( colc_row_coef.begin(), colc_row_coef.end(), NA_REAL ) ;
    NumericMatrix colc_cov_coef(1,1);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;
    NumericMatrix rowc_colc_coef(1,1);
    std::fill( rowc_colc_coef.begin(), rowc_colc_coef.end(), NA_REAL ) ;
    int CG = 0;

    rcpp_unpack(model_num, invect, param_lengths, mu, phi, rowc_coef, colc_coef,
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

    // Adjust very small values of pi_v to avoid errors when calculating log(pi_v)
    for (rr=0; rr < RG; rr++) {
        if (pi_v[rr] < epsilon) pi_v[rr] = epsilon;
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
                if (is_true(all(is_finite(yval))) && is_true(all(!is_nan(yval)))) {
                    linear_part = rcpp_linear_part(ydf, rowc_mm, colc_mm, cov_mm,
                                                   param_lengths,
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

                    theta = rcpp_theta_from_linear(model_num, linear_part, ymatij_idx, mu, phi, q, epsilon);

                    log_thetaymat = log(theta);
                } else {
                    log_thetaymat = 0;
                }
                // Rcout << "The value of log_thetaymat : " << log_thetaymat << "\n";

                // Rcout << "The value of llc component : " << ppr_m(ii,rr)*log_thetaymat << "\n";
                llc += ppr_m(ii,rr)*log_thetaymat;
            }
        }

        if (!partial) {
            for (ii=0; ii < n; ii++) {
                for (rr=0; rr < RG; rr++) {
                    llc += ppr_m(ii,rr)*log(pi_v[rr]);
                }
            }
        }

        logl = llc;
    } else {
        logl = 0;

        NumericVector log_components (RG);
        double log_sumoverR = 0;
        double theta;
        double log_theta;

        for (ii=0; ii < n; ++ii) {
            // Rcout << "The R-based value of ii : " << ii+1 << "\n";

            log_components.fill(0);

            for (rr=0; rr < RG; rr++) {
                log_components[rr] = log(pi_v[rr]);

                for (ij=0; ij < ydf.nrow(); ij++) {
                    if (ydf(ij,1) == ii+1) {

                        jj = ydf(ij,2)-1;
                        // Rcout << "The R-based value of jj : " << jj+1 << "\n";
                        yval = ydf(ij,0);
                        if (is_true(all(is_finite(yval))) && is_true(all(!is_nan(yval)))) {
                            linear_part = rcpp_linear_part(ydf, rowc_mm, colc_mm, cov_mm,
                                                           param_lengths,
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

                            theta = rcpp_theta_from_linear(model_num, linear_part, ymatij_idx, mu, phi, q, epsilon);
                            log_theta = log(theta);
                            if (!NumericVector::is_na(log_theta) &&
                                !std::isnan(log_theta) &&
                                std::isfinite(log_theta)) {
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
                        const int & model_num,
                        const NumericMatrix & ydf,
                        const NumericMatrix & rowc_mm,
                        const NumericMatrix & colc_mm,
                        const NumericMatrix & cov_mm,
                        const NumericMatrix & ppr_m,
                        const NumericMatrix & ppc_m,
                        NumericVector pi_v,
                        NumericVector kappa_v,
                        const IntegerVector & param_lengths,
                        const int & RG, const int & CG,
                        const int & p, const int & n, const int & q,
                        const double & epsilon,
                        const bool & constraint_sum_zero, const bool & partial,
                        const bool & incomplete, double & llc) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: param_lengths["rowc"], else you will get a
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
    if (param_lengths["rowc_cov"] > 0) {
        nrowccov = param_lengths["rowc_cov"]/RG;
    }
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;

    int ncolccov = 0;
    if (param_lengths["colc_cov"] > 0) {
        ncolccov = param_lengths["colc_cov"]/CG;
    }
    // Rcout << "The value of ncolccov: " << ncolccov << "\n";
    NumericMatrix colc_cov_coef(CG,ncolccov);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;

    NumericVector cov_coef(param_lengths["cov"],NA_REAL);

    rcpp_unpack(model_num, invect, param_lengths, mu, phi, rowc_coef, colc_coef,
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

    double pi_kappa_component;

    // Adjust very small values of pi_v to avoid errors when calculating log(pi_v)
    for (rr=0; rr < RG; rr++) {
        if (pi_v[rr] < epsilon) pi_v[rr] = epsilon;
    }
    // Similarly for kappa_v
    for (cc=0; cc < CG; cc++) {
        if (kappa_v[cc] < epsilon) kappa_v[cc] = epsilon;
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
                    if (is_true(all(is_finite(yval))) && is_true(all(!is_nan(yval)))) {
                        ymatij_idx = ydf(ij,0)-1;

                        // Rcout << "The value of ymatij_idx : " << ymatij_idx << "\n";

                        linear_part = rcpp_linear_part(ydf, rowc_mm, colc_mm, cov_mm,
                                                       param_lengths,
                                                       rowc_coef, colc_coef,
                                                       rowc_colc_coef,
                                                       row_coef, col_coef,
                                                       rowc_col_coef, colc_row_coef,
                                                       rowc_cov_coef, colc_cov_coef,
                                                       cov_coef,
                                                       RG, CG, p, n, q,
                                                       nrowccov, ncolccov,
                                                       rr, cc, ij, ii, jj);

                        theta = rcpp_theta_from_linear(model_num, linear_part, ymatij_idx, mu, phi, q, epsilon);
                        // Rcout << "The value of theta : " << theta << " " << temp << "\n";

                        log_thetaymat = log(theta);
                    } else {
                        log_thetaymat = 0;
                    }
                    // Rcout << "The value of log_thetaymat : " << log_thetaymat << "\n";
                    // Rcout << "The value of llc component : " << ppr_m(ii,rr)*log_thetaymat*ppc_m(jj,cc) << "\n";
                    llc += ppr_m(ii,rr)*log_thetaymat*ppc_m(jj,cc);
                }
            }
        }

        if (!partial) {
            // Rcout << "The value of pi_v : " << pi_v << "\n";
            // Rcout << "The value of kappa_v : " << kappa_v << "\n";

            // Rcout << "The value of ppr : " << ppr_m << "\n";
            // Rcout << "The value of ppc : " << ppc_m << "\n";

            for (ii=0; ii < n; ii++) {
                for (rr=0; rr < RG; rr++) {
                    pi_kappa_component = ppr_m(ii,rr)*log(pi_v[rr]);
                    if (!std::isnan(pi_kappa_component)) {
                        llc += pi_kappa_component;
                    }
                }
            }
            // Rcout << "The value of llc : " << llc << "\n";
            for (jj=0; jj < p; jj++) {
                for (cc=0; cc < CG; cc++) {
                    pi_kappa_component = ppc_m(jj,cc)*log(kappa_v[cc]);
                    if (!std::isnan(pi_kappa_component)) {
                        llc += pi_kappa_component;
                    }
                }
            }
            // Rcout << "The value of llc : " << llc << "\n";
        }

        logl = llc;
    }

    if (incomplete) {
        double logl_correction = 0;

        for (ii=0; ii < n; ii++) {
            for (rr=0; rr < RG; rr++) {
                if (ppr_m(ii,rr) > 0) {
                    // Note that the "n*" part is the sum of xhat_jc over j and c,
                    // because the original term is sum_i sum_j sum_r sum_c zhat_ir xhat_jc log(zhat_ir),
                    // but the xhat_jc bit simplifies to p
                    logl_correction += p*ppr_m(ii,rr)*log(ppr_m(ii,rr));
                }
            }
        }

        for (jj=0; jj < p; jj++) {
            for (cc=0; cc < CG; cc++) {
                if (ppc_m(jj,cc) > 0) {
                    // Note that the "p*" part is the sum of zhat_ir over i and r
                    // because the original term is sum_i sum_j sum_r sum_c zhat_ir xhat_jc log(xhat_jc)
                    // but the zhat_ir bit simplifies to n
                    logl_correction += n*ppc_m(jj,cc)*log(ppc_m(jj,cc));
                }
            }
        }

        if (!NumericVector::is_na(logl_correction) &&
            !std::isnan(logl_correction) &&
            std::isfinite(logl_correction)) {
            logl = llc - logl_correction;
        } else {
            logl = R_NegInf;
        }
        // Rcout << "The value of logl : " << logl << "\n";
    }

    return logl;
}

// [[Rcpp::export]]
NumericMatrix rcpp_Rcluster_Estep(const NumericVector & invect,
                                  const int & model_num,
                                  const NumericMatrix & ydf,
                                  const NumericMatrix & rowc_mm,
                                  const NumericMatrix & colc_mm,
                                  const NumericMatrix & cov_mm,
                                  NumericVector pi_v,
                                  const IntegerVector & param_lengths,
                                  const int & RG, const int & p, const int & n, const int & q,
                                  const double & epsilon,
                                  const bool & constraint_sum_zero) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: param_lengths["rowc"], else you will get a
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

    int nrowccov = param_lengths["rowc_cov"]/RG;
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;
    int ncolccov = 0;

    NumericVector cov_coef(param_lengths["cov"],NA_REAL);

    NumericVector colc_coef(1,NA_REAL);
    NumericVector row_coef(1,NA_REAL);
    NumericMatrix colc_row_coef(1,1);
    std::fill( colc_row_coef.begin(), colc_row_coef.end(), NA_REAL ) ;
    NumericMatrix colc_cov_coef(1,1);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;
    NumericMatrix rowc_colc_coef(1,1);
    std::fill( rowc_colc_coef.begin(), rowc_colc_coef.end(), NA_REAL ) ;
    int CG = 0;

    rcpp_unpack(model_num, invect, param_lengths, mu, phi, rowc_coef, colc_coef,
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

    // Adjust very small values of pi_v to avoid errors when calculating log(pi_v)
    for (rr=0; rr < RG; rr++) {
        if (pi_v[rr] < epsilon) pi_v[rr] = epsilon;
    }

    for (ij=0; ij < ydf.nrow(); ij++) {

        ii = ydf(ij,1)-1;
        jj = ydf(ij,2)-1;

        for (rr=0; rr < RG; rr++) {
            if (ppm_started(ii,rr) == 0) {
                ppm_raw(ii,rr) = log(pi_v[rr]);
                ppm_started(ii,rr) = 1;
            }

            yval = ydf(ij,0);
            if (is_true(all(is_finite(yval))) && is_true(all(!is_nan(yval)))) {
                linear_part = rcpp_linear_part(ydf, rowc_mm, colc_mm, cov_mm,
                                               param_lengths,
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

                theta = rcpp_theta_from_linear(model_num, linear_part, ymatij_idx, mu, phi, q, epsilon);

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
                                   const int & model_num,
                                   const NumericMatrix & ydf,
                                   const NumericMatrix & rowc_mm,
                                   const NumericMatrix & colc_mm,
                                   const NumericMatrix & cov_mm,
                                   NumericVector pi_v,
                                   NumericVector kappa_v,
                                   const IntegerVector & param_lengths,
                                   const int & RG, const int & CG,
                                   const int & p, const int & n, const int & q,
                                   const double & epsilon,
                                   const bool & constraint_sum_zero,
                                   const bool & row_clusters) {

    // BE CAREFUL! You CANNOT use single quotes to fetch named element of a
    // vector, only double quotes: param_lengths["rowc"], else you will get a
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
    if (param_lengths["rowc_cov"] > 0) {
        nrowccov = param_lengths["rowc_cov"]/RG;
    }
    NumericMatrix rowc_cov_coef(RG,nrowccov);
    std::fill( rowc_cov_coef.begin(), rowc_cov_coef.end(), NA_REAL ) ;

    int ncolccov = 0;
    if (param_lengths["colc_cov"] > 0) {
        ncolccov = param_lengths["colc_cov"]/CG;
    }
    // Rcout << "The value of ncolccov: " << ncolccov << "\n";
    NumericMatrix colc_cov_coef(CG,ncolccov);
    std::fill( colc_cov_coef.begin(), colc_cov_coef.end(), NA_REAL ) ;

    NumericVector cov_coef(param_lengths["cov"],NA_REAL);

    rcpp_unpack(model_num, invect, param_lengths, mu, phi, rowc_coef, colc_coef,
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

    // Adjust very small values of pi_v to avoid errors when calculating log(pi_v)
    for (rr=0; rr < RG; rr++) {
        if (pi_v[rr] < epsilon) pi_v[rr] = epsilon;
    }
    // Similarly for kappa_v
    for (cc=0; cc < CG; cc++) {
        if (kappa_v[cc] < epsilon) kappa_v[cc] = epsilon;
    }

    for (ij=0; ij < ydf.nrow(); ij++) {

        ii = ydf(ij,1)-1;
        jj = ydf(ij,2)-1;
        // Rcout << "The R-based value of ii : " << ii+1 << "\n";
        // Rcout << "The R-based value of jj : " << jj+1 << "\n";

        if (row_clusters) {

            for (rr=0; rr < RG; rr++) {
                // Need this step because we are looping over ij, but for each
                // ii and rr we only want to include one log(pi_v[rr]), instead of
                // adding one for each jj
                if (ppm_started(ii,rr) == 0) {
                    ppm_raw(ii,rr) = log(pi_v[rr]);
                    ppm_started(ii,rr) = 1;
                }

                yval = ydf(ij,0);
                if (is_true(all(is_finite(yval))) && is_true(all(!is_nan(yval)))) {
                    ymatij_idx = ydf(ij,0)-1;
                    sum_pikappa_theta = 0;

                    for (cc=0; cc < CG; cc++) {
                        linear_part = rcpp_linear_part(ydf, rowc_mm, colc_mm, cov_mm,
                                                       param_lengths,
                                                       rowc_coef, colc_coef,
                                                       rowc_colc_coef,
                                                       row_coef, col_coef,
                                                       rowc_col_coef, colc_row_coef,
                                                       rowc_cov_coef, colc_cov_coef,
                                                       cov_coef,
                                                       RG, CG, p, n, q,
                                                       nrowccov, ncolccov,
                                                       rr, cc, ij, ii, jj);

                        theta = rcpp_theta_from_linear(model_num, linear_part, ymatij_idx, mu, phi, q, epsilon);

                        sum_pikappa_theta += kappa_v[cc]*theta;
                    }

                    // Rcout << "The value of the additional log term : " << log(sum_pikappa_theta) << "\n";
                    ppm_raw(ii,rr) += log(sum_pikappa_theta);
                }
            }
        } else {
            for (cc=0; cc < CG; cc++) {
                // Need this step because we are looping over ij, but for each
                // jj and cc we only want to include one log(kappa_v[cc]),
                // instead of adding one for each ii
                if (ppm_started(jj,cc) == 0) {
                    ppm_raw(jj,cc) = log(kappa_v[cc]);
                    ppm_started(jj,cc) = 1;
                }

                yval = ydf(ij,0);
                if (is_true(all(is_finite(yval))) && is_true(all(!is_nan(yval)))) {
                    ymatij_idx = ydf(ij,0)-1;
                    sum_pikappa_theta = 0;

                    for (rr=0; rr < RG; rr++) {
                        linear_part = rcpp_linear_part(ydf, rowc_mm, colc_mm, cov_mm,
                                                       param_lengths,
                                                       rowc_coef, colc_coef,
                                                       rowc_colc_coef,
                                                       row_coef, col_coef,
                                                       rowc_col_coef, colc_row_coef,
                                                       rowc_cov_coef, colc_cov_coef,
                                                       cov_coef,
                                                       RG, CG, p, n, q,
                                                       nrowccov, ncolccov,
                                                       rr, cc, ij, ii, jj);

                        theta = rcpp_theta_from_linear(model_num, linear_part, ymatij_idx, mu, phi, q, epsilon);

                        sum_pikappa_theta += pi_v[rr]*theta;
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


