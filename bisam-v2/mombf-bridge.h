//
// Created by jakob on 3/24/25.
//

#ifndef MOMBF_BRIDGE_H
#define MOMBF_BRIDGE_H

#include <complex.h>
#include "lib/lasso/LassoRegression.h"

#include "mombf/mombf/src/modelSel_regression.h"

namespace MombfBridge {
    /* -------------------------------------------------------------------------- */
    /*                           ModelSelectionGibbsCI                            */
    /* -------------------------------------------------------------------------- */
    // The following code is a modified version of the original modelSelectionGibbsCI. The original function
    // prepares the R (SEXP) objects for use with cpp. Since we already have the cpp objects ready not much
    // is needed here. We mainly make sure all our parameters are in a form to call modelSelectionGibbs.

    void modelSelectionGibbsCI(const arma::vec &SpostModeini,
                               double SpostModeiniProb,
                               int Sknownphi,
                               int Sfamily,
                               int SpriorCoef,
                               int SpriorGroup,
                               int Sniter,
                               int Sthinning,
                               int Sburnin,
                               int Sndeltaini,
                               arma::Col<int> &Sdeltaini,
                               arma::vec &Sincludevars,
                               int Sn,
                               int Sp,
                               arma::vec &Sy,
                               int Suncens,
                               double Ssumy2,
                               double Ssumy,
                               double Ssumlogyfact,
                               arma::mat &Sx,
                               arma::vec &Scolsumsx,
                               bool ShasXtX,
                               arma::mat &SXtX,
                               arma::rowvec &SytX,
                               int Smethod,
                               int Sadjoverdisp,
                               int Shesstype,
                               int SoptimMethod,
                               int Soptim_maxit,
                               arma::vec Sthinit,
                               int Susethinit,
                               int SB,
                               double Salpha,
                               double Slambda,
                               double Sphi,
                               double Stau,
                               double Staugroup,
                               double Staualpha,
                               int Sfixatanhalpha,
                               int Sr,
                               int SpriorDelta,
                               double SprDeltap,
                               int SparprDeltap,
                               int SpriorConstr,
                               double SprConstrp,
                               int SparprConstrp,
                               int *Sgroups,
                               int Sngroups,
                               arma::uvec &Snvaringroup,
                               arma::ivec &Sconstraints,
                               arma::ivec &Sinvconstraints,
                               int Sverbose) {
        // bool hasXtX = LOGICAL(ShasXtX)[0];

        int i, j, idxj, logscale = 1, mcmc2save, *postSample, *postMode, mycols, mycols2, *nconstraints, *
                ninvconstraints, nuncens, ngroupsconstr = 0, *isgroup, usethinit = Susethinit, priorcode;

        double offset = 0, *margpp, *postModeProb, *postProb, *ytXuncens = NULL, *thinit;
        intptrvec constraints, invconstraints;
        crossprodmat *XtX, *XtXuncens = NULL;
        struct marginalPars pars;

        // TODO Answer struct does not yet exist
        // SEXP ans;
        // PROTECT(ans = Rf_allocVector(VECSXP, 5));

        mcmc2save = floor((Sniter - Sburnin + .0) / (Sthinning + .0));
        if (Sfamily != 0) {
            mycols = mycols2 = Sp;
        } else {
            mycols  = 2 + Sp;
            mycols2 = mycols + 2;
        }

        thinit = dvector(0, mycols2 + 1);
        if (usethinit != 3) {
            for (j = 0; j <= mycols2 + 1; j++) { thinit[j] = 0; }
        } else {
            for (j = 0; j <= Sp; j++) { thinit[j] = Sthinit[j]; }
        }

        // SET_VECTOR_ELT(ans, 0, Rf_allocVector(INTSXP, mcmc2save * mycols));
        // postSample = INTEGER(VECTOR_ELT(ans, 0));
        for (j = 0; j < (mcmc2save * mycols); j++) postSample[j] = 0;

        // SET_VECTOR_ELT(ans, 1, Rf_allocVector(REALSXP, mycols2));
        // margpp = REAL(VECTOR_ELT(ans, 1));
        //
        // SET_VECTOR_ELT(ans, 2, Rf_allocVector(INTSXP, mycols));
        // postMode = INTEGER(VECTOR_ELT(ans, 2));
        for (j = 0; j < mycols; j++) { postMode[j] = SpostModeini[j]; }

        // SET_VECTOR_ELT(ans, 3, Rf_allocVector(REALSXP, 1));
        // postModeProb    = REAL(VECTOR_ELT(ans, 3));
        postModeProb[0] = SpostModeiniProb;

        // SET_VECTOR_ELT(ans, 4, Rf_allocVector(REALSXP, mcmc2save));
        // postProb = REAL(VECTOR_ELT(ans, 4));
        //
        isgroup         = ivector(0, Sp);
        nconstraints    = ivector(0, Sngroups);
        ninvconstraints = ivector(0, Sngroups);
        // TODO this is quite annoying since it uses SEXP objects which we do not have! Bridge needed
        // countConstraints(nconstraints,
        //                  &constraints,
        //                  ninvconstraints,
        //                  &invconstraints,
        //                  &ngroupsconstr,
        //                  isgroup,
        //                  &Sngroups,
        //                  (int *) Snvaringroup.memptr(),
        //                  &Sconstraints,
        //                  Sinvconstraints);

        if (ShasXtX) {
            XtX = new crossprodmat(SXtX.memptr(), Sn, Sp, true);
        } else {
            XtX = new crossprodmat(Sx.memptr(), Sn, Sp, false);
        }

        if (Suncens > 0) {
            //if there's censoring, also store t(x) %*% x and t(x) %*% y computed over uncensored observations
            int n       = Sn, *uncens       = &Suncens; // JAKOB added & in front of Suncens since its already int.
            double *pty = Sy.memptr(), *ptx = Sx.memptr();
            for (nuncens = 0; (nuncens < n) && (uncens[nuncens] == 1); nuncens++) {
            }
            XtXuncens = new crossprodmat(Sx.memptr(), Sn, Sp, false, nuncens, 0);
            ytXuncens = dvector(0, Sp);
            for (j = 0; j < Sp; j++) {
                for (i = 0, ytXuncens[j] = 0, idxj = j * n; i < nuncens; i++) {
                    ytXuncens[j] += pty[i] * ptx[i + idxj];
                }
            }
        } else { nuncens = Sn; }

        set_marginalPars(&pars,
                         &Sfamily,
                         &Sn,
                         &nuncens,
                         &Sp,
                         Sy.memptr(),
                         &Suncens,
                         &Ssumy2,
                         &Ssumy,
                         &Ssumlogyfact,
                         Sx.memptr(),
                         Scolsumsx.memptr(),
                         XtX,
                         SytX.memptr(),
                         &Smethod,
                         &Sadjoverdisp,
                         &Shesstype,
                         &SoptimMethod,
                         &Soptim_maxit,
                         &usethinit,
                         thinit,
                         &SB,
                         &Salpha,
                         &Slambda,
                         &Sknownphi,
                         &Sphi,
                         &Stau,
                         &Staugroup,
                         &Staualpha,
                         (double *) Sfixatanhalpha,
                         &Sr,
                         &SprDeltap,
                         (double *) SparprDeltap,
                         &SprConstrp,
                         (double *) SparprConstrp,
                         &logscale,
                         &offset,
                         Sgroups,
                         isgroup,
                         &Sngroups,
                         &ngroupsconstr,
                         // TODO not super clean.
                         (int *) Snvaringroup.memptr(),
                         nconstraints,
                         ninvconstraints,
                         XtXuncens,
                         ytXuncens);

        priorcode      = mspriorCode(&SpriorCoef, &SpriorGroup, &pars);
        pars.priorcode = &priorcode;

        modelSelectionGibbs(postSample,
                            margpp,
                            postMode,
                            postModeProb,
                            postProb,
                            &SpriorDelta,
                            &SpriorConstr,
                            &Sniter,
                            &Sthinning,
                            &Sburnin,
                            &Sndeltaini,
                            Sdeltaini.memptr(),
                            (int *) Sincludevars.memptr(),
                            &constraints,
                            &invconstraints,
                            &Sverbose,
                            &pars);

        free_dvector(thinit, 0, mycols2 + 1);
        free_ivector(isgroup, 0, Sp);
        free_ivector(nconstraints, 0, Sngroups);
        free_ivector(ninvconstraints, 0, Sngroups);
        delete XtX;
        // return ans;
    }

    /* -------------------------------------------------------------------------- */
    /*                                 getthinit                                  */
    /* -------------------------------------------------------------------------- */
    // The following code provides 2 functions and an enum to emulate the getthinit and initParameters
    // R functions.


    // Enum to replace the strings used in the original R code
    enum InitparType {
        None,
        MLE,
        MLE_aisg,
        L1,
        L2_aisgd,
        Auto,
    };

    // The original function used the initpar as a string (always a string and never vec here!). Instead
    // we use our enum. For now we assume family to always be 1 (=normal!). Thus some code is not translated
    arma::vec initParameters(const arma::vec &y, const arma::mat &x, int family,
                             InitparType initpar) {
        int n = y.n_elem;
        int p = x.n_cols;

        if (initpar == InitparType::Auto) {
            if (p <= n / 2) {
                initpar = InitparType::MLE;
            } else {
                initpar = InitparType::L1;
            }
        }

        // I think these lines should never be needed. Setting the family here should not do anything since we
        // are hard coding it later. And initpar should never be none here. 
        // if (!(family % in % c('binomial', 'poisson'))) family = 'gaussian'
        // if (initpar == 'none') {
        //     ans = rep(0, p)
        // }

        // TODO if needed the lasso code can probably easily be adapted to accept arma types directly
        // Convert Armadillo data to the format expected by LassoRegression
        std::vector<std::vector<double> > samples(n, std::vector<double>(p));
        std::vector<double> target(n);

        for (int i = 0; i < n; ++i) {
            target[i] = y(i);
            for (int j = 0; j < p; ++j) {
                samples[i][j] = x(i, j);
            }
        }

        // Furthermore I think that we ALWAYS hit the L1 case? So we will simply use a lasso here not matter what
        // Create LassoRegression object
        LassoRegression lasso(samples, target);

        // Set tolerance and alpha for coordinate descent
        double tolerance = 0.001;
        double alpha     = 0.01; // This corresponds to the regularization strength

        // Run cyclical coordinate descent to get the weights
        double *weights = lasso.cyclicalCoordinateDescent(tolerance, alpha);

        // Convert the result (double*) to arma::vec
        arma::vec arma_weights(weights, p, false, true);

        // Clean up the dynamically allocated weights array
        delete[] weights;

        return arma_weights;
    }


    // Translation of the R function found in initParameters.R
    // Since mixing types in C++ is not really a thing and the original function uses the initpar parameter
    // as either a string or a vector, instead here initpar is only a vec. If it is NOT empty, the the value
    // is used directly. If it is empty, the type in initpar_type shall be used to determin what value
    // usethinit should be set to
    std::pair<arma::vec, int> getthinit(const arma::vec &y, const arma::mat &x, int family,
                                        const arma::vec &initpar, bool enumerate, InitparType initpar_type) {
        arma::vec thinit;
        int usethinit;

        if (!initpar.empty()) {
            thinit    = initpar;
            usethinit = 1;
            if (thinit.n_elem != x.n_cols) {
                throw std::runtime_error(
                    "x has " + std::to_string(x.n_cols) + " columns but initpar has " + std::to_string(thinit.n_elem) +
                    " elements");
            }
        } else {
            if (initpar_type == InitparType::None) {
                if (enumerate) {
                    usethinit = 0;
                } else {
                    usethinit = 1;
                }
            } else {
                usethinit = 3;
            }
            thinit = initParameters(y, x, family, initpar_type);
            // Passing an empty string as initpar to match the original logic for missing/auto
        }

        return {thinit, usethinit};
    }

    /* -------------------------------------------------------------------------- */
    /*                               modelSelection                               */
    /* -------------------------------------------------------------------------- */
    // This is a (very much simplified) version of the R function. Many parameters are hard coded since their values
    // are expected to remain the way they are for now in order for us not having to translate every single function
    // used within this one.

    void modelSelection(const arma::vec &y, const arma::mat &x, int niter, int thinning, int burnin,
                        arma::vec &deltaini_input, bool center, bool scale,
                        bool XtXprecomp, double phi, double tau, double priorSkew,
                        double prDeltap, arma::vec thinit, InitparType initpar_type) {
        int p = x.n_cols;
        int n = y.n_elem;

        /* -------------------------------- postMode -------------------------------- */
        // Since familyint == 1 for now -> postMode is a vector of length p filled with zeros
        arma::vec postMode = arma::vec().zeros(p);

        /* ------------------------------ postModeProb ------------------------------ */
        double postModeProb = 0.0; // Initialized to 0.0

        /* -------------------------------- knownphi -------------------------------- */
        // We give phi as parameter -> phi is known
        // TODO potentially change to bool
        int knownphi = 1;

        /* ------------------------------- familyint -------------------------------- */
        // Set to 1 since we have family == "normal" and !issurvival (in formatFamily())
        int familyint = 1;

        /* --------------------------------- prior ---------------------------------- */
        // Set to 1. Defined in formatmsPriorsMarg. Equal to prior = piMOM
        int prior = 1;

        /* -------------------------------- priorgr --------------------------------- */
        // Set to 1. Defined in formatmsPriorsMarg. Equal to prior = piMOM
        int priorgr = 1;

        /* -------------------------- ndeltaini & deltaini -------------------------- */
        // Translation of the respective R code:
        // ndeltaini <- as.integer(sum(deltaini | includevars))
        // deltaini <- as.integer(which(deltaini | includevars) - 1)
        arma::vec includevars = arma::vec().zeros(p);

        int ndeltaini           = arma::sum(deltaini_input || includevars);
        arma::uvec indices      = arma::find(deltaini_input || includevars);
        arma::Col<int> deltaini = arma::conv_to<arma::Col<int> >::from(indices);

        /* ------------------- ystd, sumy2, sumy, xstd, colsumsx -------------------- */
        // The following calculations depend on one another and an actual translation is needed.
        // They have all been put in this one section since they are related.
        // TODO This part may include that is not of interest if scaling is never of interest?

        // Pre define variables so they are available after the if blocks aswell
        double my;
        double sy;

        // mx <- colMeans(x)
        arma::rowvec mx = arma::mean(x, 0);

        // sx <- sqrt(colMeans(x^2) - mx^2) * sqrt(n / (n - 1))
        // TODO Theoretically the case where n = 1 would need to be caught and handled to avoid div by zero
        arma::rowvec sx = (arma::sqrt(arma::mean(arma::pow(x, 2), 0) - arma::pow(mx, 2)))
                          * std::sqrt(static_cast<double>(n) / (n - 1));

        // ct <- (sx == 0)
        // TODO very much not sure if this is correct?
        arma::uvec ct = arma::conv_to<arma::uvec>::from(sx == 0);

        // if (!center) {
        //      my <- 0
        //      mx <- rep(0, p)
        // } else {
        //      my <- mean(y)
        // }
        if (!center) {
            double my = 0.0;
            arma::vec mx(p, arma::fill::zeros);
        } else {
            double my = arma::mean(y);
        }

        // if (!scale) {
        //      sy <- 1
        //      sx <- rep(1,p)
        // } else {
        //      sy <- sd(y)
        // }
        if (!scale) {
            double sy = 1.0;
            arma::vec sx(p, arma::fill::ones);
        } else {
            double sy = arma::stddev(y);
        }

        // Since typeofvar is always "numeric" in our case so far, this code block is skipped.
        // mx[typeofvar == "factor"] <- 0
        // sx[typeofvar == "factor"] <- 1

        // Since for now always false, this section of code is skipped
        // if (!(outcometype %in% c("Continuous", "Survival"))) {
        //   my <- 0
        //   sy <- 1
        // }

        // ystd <- (y - my) / sy
        arma::vec ystd = (y - my) / sy;

        // xstd <- x
        arma::mat xstd = x;

        // xstd[, !ct] <- t((t(x[, !ct]) - mx[!ct]) / sx[!ct])
        // In R, this standardizes columns with non-zero standard deviation
        // In Armadillo, we'll loop through columns and only standardize those with ct(j) == 0
        for (size_t j = 0; j < p; ++j) {
            if (ct(j) == 0) {
                // If sx(j) is not 0 (equivalent to !ct in R)
                xstd.col(j) = (x.col(j) - mx(j)) / sx(j);
            }
        }

        // sumy2 <- as.double(sum(ystd^2))
        // Using dot product which is more efficient for sum of squares
        double sumy2 = arma::dot(ystd, ystd);

        // sumy <- as.double(sum(ystd))
        double sumy = arma::sum(ystd);

        // ytX <- as.vector(matrix(ystd, nrow = 1) %*% xstd)
        // In Armadillo, we can multiply a column vector's transpose with a matrix to get row vector
        arma::rowvec ytX = ystd.t() * xstd;

        // colsumsx <- as.double(colSums(xstd))
        // 0 indicates summing along rows to get column sums
        arma::vec colsumsx = arma::sum(xstd, 0).t();

        // if (XtXprecomp) { XtX <- t(xstd) %*% xstd; hasXtX <- as.logical(TRUE) }
        // else { XtX <- double(0); hasXtX <- as.logical(FALSE) }
        arma::mat XtX;
        bool hasXtX;
        if (XtXprecomp) {
            // Compute the cross-product matrix
            XtX    = xstd.t() * xstd;
            hasXtX = true;
        } else {
            // Create an empty matrix instead of a single double
            XtX.reset();
            hasXtX = false;
        }

        /* ------------------------------ sumlogyfact ------------------------------- */
        // Hard coded since family is not 22 sumlogyfact <- as.double(0)
        double sumlogyfact = 0;

        /* --------------------------------- uncens --------------------------------- */
        // Hard coded to 0 since "formula" %in% class(y) = False and "Surv" %in% class(y) in
        // formatInputData
        int uncens = 0;

        /* --------------------------------- method --------------------------------- */
        // Fixed to int 0 in formats method. Because opti method = auto and outcome type = glm
        int method = 0;

        /* ------------------------------ adj_overdisp ------------------------------ */
        // Fixed to 1 because of formats method. adj.overdisp <- as.integer(ifelse(adj.overdisp == "none", 0,
        // ifelse(adj.overdisp == "intercept", 1, 2))). We set adj.overdisp='intercept'.
        int adj_overdisp = 1;

        /* -------------------------------- hesstype -------------------------------- */
        // Fixed to 1 because of formats method. hesstype <- as.integer(ifelse(hess == "asympDiagAdj", 2, 1))
        int hesstype = 1;

        /* ------------------------------ optimMethod ------------------------------- */
        // Fixed to 2 because of formats method. We set optimMethod="auto".
        int optimMethod = 2;

        /* ------------------------------ optim_maxit ------------------------------- */
        // Default is set to 10? Not really a default but we set it in our list of "default" parameters?
        int optim_maxit = 10;

        /* -------------------------- thinit and usethinit -------------------------- */
        // Thinit uses a lasso for initialization - put into its own function for clearity
        int usethinit;

        auto res  = getthinit(y, x, familyint, thinit, false, initpar_type);
        thinit    = res.first;
        usethinit = res.second;


        /* ----------------------------------- B ------------------------------------ */
        // Default parameter set to 100000
        int B = 100000;

        /* --------------------------------- alpha ---------------------------------- */
        // Fixed to 0.01 - parameter set by us. alpha <- as.double(priorVar@priorPars["alpha"])
        double alpha = 0.01;

        /* --------------------------------- lambda --------------------------------- */
        // Fixed to 0.01 - parameter set by us. lambda <- as.double(priorVar@priorPars["lambda"])
        double lambda = 0.01;

        /* ---------------------------------- phi ----------------------------------- */
        // Set to be s2_i by us. Parameter of this function no further calculations.

        /* ---------------------------------- tau ----------------------------------- */
        // priorCoef = imomprior(tau = tau); originally. We just add a function paramter to get tau

        /* -------------------------------- taugroup -------------------------------- */
        // We have missing(priorgroup) == True ---> priorGroup <- priorCoef
        // Thus taugroup == tau? Hastaugroupstd == false ---> taugroup <- as.double(priorGroup@priorPars["tau"])
        double taugroup = tau;

        /* -------------------------------- taualpha -------------------------------- */
        // Since "msPriorSpec" %in% class(priorSkew) == true ---> taualpha <- as.double(priorSkew@priorPars["tau"])
        double taualpha = priorSkew;

        /* ----------------------------- fixatanhalpha ------------------------------ */
        // Since "msPriorSpec" %in% class(priorSkew) == true ---> fixatanhalpha <- as.double(-10000)
        int fixatanhalpha = -10000;

        /* ----------------------------------- r ------------------------------------ */
        // Fixed to 1. In formatmsPriorsMarg for piMom
        int r = 1;

        /* -------------------------------- prDelta --------------------------------- */
        int prDelta = 1;

        /* -------------------------------- prDeltap -------------------------------- */
        // Parameter priorDelta@priorPars[["p"]]

        /* ------------------------------ parprDeltap ------------------------------- */
        // parprDeltap <- as.double(length(prDeltap))
        int parprDeltap = 1;

        /* -------------------------------- prConstr -------------------------------- */
        // Fixed to one in formatmsPriorsMarg
        int prConstr = 1;

        /* ------------------------------- prConstrp -------------------------------- */
        // priorConstraints < -defaultpriorConstraints(priorDelta)
        // prConstrp <- as.double(priorConstraints@priorPars[["p"]])
        // Fixed to 0.5 in our case
        double prConstrp = 0.5;

        /* ------------------------------ parprConstrp ------------------------------ */
        // parprConstrp <- as.double(length(prConstrp))
        int parprConstrp = 1;

        /* --------------------------------- groups --------------------------------- */
        // Vector with numbers 0 to 20. Groups is modiefied (compared to default parameter) in
        // codeGroupsAndConstraints
        arma::uvec groups = arma::regspace<arma::uvec>(0, p - 1);

        /* -------------------------------- ngroups --------------------------------- */
        // = p
        int ngroups = p;

        /* ------------------------------ nvaringroup ------------------------------- */
        // nvaringroup <- as.integer(rep(1, p))
        arma::uvec nvaringroup = arma::regspace<arma::uvec>(1, p);

        /* ------------------------------ constraints ------------------------------- */
        // I am not 100% sure about constraints and invconstraints. Further checking may be needed
        arma::ivec constraints = arma::ivec().zeros(p);

        /* ----------------------------- invconstraints ----------------------------- */
        // I am not 100% sure about constraints and invconstraints. Further checking may be needed
        arma::ivec invconstraints = arma::ivec().zeros(p);

        /* -------------------------------- verbose --------------------------------- */
        int verbose = 1;

        modelSelectionGibbsCI(postMode,
                              postModeProb,
                              knownphi,
                              familyint,
                              prior,
                              priorgr,
                              niter,
                              thinning,
                              burnin,
                              ndeltaini,
                              deltaini,
                              includevars,
                              n,
                              p,
                              ystd,
                              uncens,
                              sumy2,
                              sumy,
                              sumlogyfact,
                              xstd,
                              colsumsx,
                              hasXtX,
                              XtX,
                              ytX,
                              method,
                              adj_overdisp,
                              hesstype,
                              optimMethod,
                              optim_maxit,
                              thinit,
                              usethinit,
                              B,
                              alpha,
                              lambda,
                              phi,
                              tau,
                              taugroup,
                              taualpha,
                              fixatanhalpha,
                              r,
                              prDelta,
                              prDeltap,
                              parprDeltap,
                              prConstr,
                              prConstrp,
                              parprConstrp,
                              (int *) groups.memptr(),
                              ngroups,
                              nvaringroup,
                              constraints,
                              invconstraints,
                              verbose);
    }
}
#endif //MOMBF_BRIDGE_H
