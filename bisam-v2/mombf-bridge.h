//
// Created by jakob on 3/24/25.
//

#ifndef MOMBF_BRIDGE_H
#define MOMBF_BRIDGE_H

struct FormattedData {
    arma::vec y;
    arma::mat x;
    std::string formula;
    arma::vec splineDegree;
    arma::vec groups;
    bool hasgroups;
    arma::vec isgroups;
    arma::vec constraints;
    std::string outcometype;
    arma::vec uncens;
    arma::vec ordery;
    std::vector<std::string> typeofvar;
};

FormattedData formatInputdata_cpp(const arma::vec &y_in, const arma::mat &x_in,
                                  const arma::mat &data, const arma::mat &smoothterms,
                                  int nknots, const std::string &family) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'formatInputdata' function's behavior.
    FormattedData formatted_data;
    formatted_data.y       = y_in;
    formatted_data.x       = x_in;
    formatted_data.formula = "placeholder_formula";
    formatted_data.splineDegree.fill(1); // Example
    return formatted_data;
}

// Placeholder for R's 'codeGroupsAndConstraints' function.
struct GroupConstraintInfo {
    int ngroups;
    arma::vec constraints;
    arma::vec invconstraints;
    arma::vec nvaringroup;
    arma::vec groups;
};

GroupConstraintInfo codeGroupsAndConstraints_cpp(int p, const arma::vec &groups_in,
                                                 const arma::vec &constraints_in) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'codeGroupsAndConstraints' function's behavior.
    GroupConstraintInfo info;
    info.ngroups     = arma::max(groups_in) + 1;
    info.constraints = constraints_in;
    info.invconstraints.zeros(constraints_in.n_elem);
    info.nvaringroup.ones(info.ngroups);
    info.groups = groups_in;
    return info;
}

// Placeholder for R's 'defaultmom' function.
struct DefaultPrior {
    arma::vec priorCoef;
    arma::vec priorVar;
};

DefaultPrior defaultmom_cpp(const std::string &outcometype, const std::string &family) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'defaultmom' function's behavior.
    DefaultPrior prior;
    prior.priorCoef.fill(0); // Example
    prior.priorVar.fill(1);  // Example
    return prior;
}

// Placeholder for R's 'groupzellnerprior' function.
arma::vec groupzellnerprior_cpp(int tau) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'groupzellnerprior' function's behavior.
    return arma::ones<arma::vec>(1); // Example
}

// Placeholder for R's 'modelbbprior' function.
arma::vec modelbbprior_cpp(double a, double b) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'modelbbprior' function's behavior.
    return arma::ones<arma::vec>(1); // Example
}

// Placeholder for R's 'igprior' function.
arma::vec igprior_cpp(double a, double b) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'igprior' function's behavior.
    return arma::ones<arma::vec>(1); // Example
}

// Placeholder for R's 'momprior' function.
arma::vec momprior_cpp(double tau) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'momprior' function's behavior.
    return arma::ones<arma::vec>(1); // Example
}

// Placeholder for R's 'getthinit' function.
struct ThinitInfo {
    bool usethinit;
    arma::vec thinit;
};

ThinitInfo getthinit_cpp(const arma::vec &y, const arma::mat &x, const std::string &family,
                         const std::string &initpar, bool enumerate) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'getthinit' function's behavior.
    ThinitInfo info;
    info.usethinit = false;
    info.thinit.zeros(1);
    return info;
}

// Placeholder for R's 'formatmsMethod' function.
struct MethodInfo {
    std::string optimMethod;
    int optim_maxit;
    std::string adj_overdisp;
    std::string hesstype;
    std::string method;
};

MethodInfo formatmsMethod_cpp(const std::string &method_in, bool usethinit,
                              const std::string &optimMethod_in, int optim_maxit_in,
                              const arma::vec &priorCoef, const arma::vec &priorGroup,
                              int knownphi, const std::string &outcometype,
                              const std::string &family, bool hasgroups,
                              const std::string &adj_overdisp_in, const std::string &hess_in) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'formatmsMethod' function's behavior.
    MethodInfo info;
    info.optimMethod  = optimMethod_in;
    info.optim_maxit  = optim_maxit_in;
    info.adj_overdisp = adj_overdisp_in;
    info.hesstype     = hess_in;
    info.method       = method_in;
    return info;
}

// Placeholder for R's 'formatFamily' function.
struct FamilyInfo {
    int familyint;
    int familygreedy;
};

FamilyInfo formatFamily_cpp(const std::string &family, bool issurvival) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'formatFamily' function's behavior.
    FamilyInfo info;
    if (family == "normal") {
        info.familyint    = 0;
        info.familygreedy = 0;
    } else {
        info.familyint    = -1; // Indicate other families
        info.familygreedy = -1;
    }
    return info;
}

// Placeholder for R's 'formatmsPriorsMarg' function.
struct MarginalPriorInfo {
    double r;
    arma::vec prior;
    arma::vec priorgr;
    double tau;
    double taugroup;
    double alpha;
    double lambda;
    double taualpha;
    bool fixatanhalpha;
    arma::vec priorCoef;
    arma::vec priorGroup;
};

MarginalPriorInfo formatmsPriorsMarg_cpp(const arma::vec &priorCoef_in,
                                         const arma::vec &priorGroup_in,
                                         const arma::vec &priorVar,
                                         const arma::vec &priorSkew, int n) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'formatmsPriorsMarg' function's behavior.
    MarginalPriorInfo info;
    info.r             = 1.0;
    info.prior         = priorCoef_in;
    info.priorgr       = priorGroup_in;
    info.tau           = 1.0;
    info.taugroup      = 1.0;
    info.alpha         = 1.0;
    info.lambda        = 1.0;
    info.taualpha      = 1.0;
    info.fixatanhalpha = false;
    info.priorCoef     = priorCoef_in;
    info.priorGroup    = priorGroup_in;
    return info;
}

// Placeholder for R's 'defaultpriorConstraints' function.
arma::vec defaultpriorConstraints_cpp(const arma::vec &priorDelta,
                                      const arma::vec &priorConstraints) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'defaultpriorConstraints' function's behavior.
    if (priorConstraints.is_empty()) {
        return priorDelta;
    } else {
        return priorConstraints;
    }
}

// Placeholder for R's 'formatmsPriorsModel' function.
struct ModelPriorInfo {
    arma::vec prDelta;
    arma::vec prDeltap;
    arma::vec parprDeltap;
    arma::vec prConstr;
    arma::vec prConstrp;
    arma::vec parprConstrp;
};

ModelPriorInfo formatmsPriorsModel_cpp(const arma::vec &priorDelta,
                                       const arma::vec &priorConstraints,
                                       const arma::vec &constraints) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'formatmsPriorsModel' function's behavior.
    ModelPriorInfo info;
    info.prDelta      = priorDelta;
    info.prDeltap     = arma::ones<arma::vec>(priorDelta.n_elem);
    info.parprDeltap  = arma::ones<arma::vec>(priorDelta.n_elem);
    info.prConstr     = priorConstraints;
    info.prConstrp    = arma::ones<arma::vec>(priorConstraints.n_elem);
    info.parprConstrp = arma::ones<arma::vec>(priorConstraints.n_elem);
    return info;
}

// Placeholder for R's 'greedyVarSelCI' function. This would require a significant
// amount of statistical implementation.
struct GreedyVarSelResult {
    arma::vec postMode;
    double postModeProb;
};

GreedyVarSelResult greedyVarSelCI_cpp(int knownphi, int familygreedy,
                                      const arma::vec &prior, const arma::vec &priorgr,
                                      int niterGreed, int ndeltaini, const arma::ivec &deltaini,
                                      const arma::ivec &includevars, int n, int p,
                                      const arma::vec &ystd, const arma::vec &uncens,
                                      double sumy2, double sumy, double sumlogyfact,
                                      const arma::mat &xstd, const arma::vec &colsumsx,
                                      bool hasXtX, const arma::mat &XtX, const arma::vec &ytX,
                                      const std::string &method, const std::string &adj_overdisp,
                                      const std::string &hesstype, const std::string &optimMethod,
                                      int optim_maxit, const arma::vec &thinit, bool usethinit,
                                      int B, double alpha, double lambda, const arma::vec &phi,
                                      double tau, double taugroup, double taualpha,
                                      bool fixatanhalpha, double r, const arma::vec &prDelta,
                                      const arma::vec &prDeltap, const arma::vec &parprDeltap,
                                      const arma::vec &prConstr, const arma::vec &prConstrp,
                                      const arma::vec &parprConstrp, const arma::vec &groups,
                                      int ngroups, const arma::vec &nvaringroup,
                                      const arma::vec &constraints, const arma::vec &invconstraints,
                                      int verbose) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'greedyVarSelCI' function's behavior and the statistical model.
    GreedyVarSelResult result;
    result.postMode.zeros(p);
    result.postModeProb = 0.0;
    return result;
}

// Placeholder for R's 'cv.ncvreg' and 'ncvreg' functions. These are related to
// penalized regression and would require external library integration or
// reimplementation.
// For now, let's assume a simplified initialization if initSearch is "SCAD".
arma::ivec initialize_scad_cpp(const arma::mat &xstd_no_ct, const arma::vec &ystd_centered,
                               int p_no_ct, const arma::ivec &includevars_no_ct, int verbose) {
    arma::ivec deltaini(p_no_ct, arma::fill::zeros);
    // Placeholder for SCAD initialization logic
    if (verbose) {
        std::cout << "Initializing via SCAD cross-validation... Done" << std::endl;
    }
    return deltaini;
}

// Placeholder for R's 'modelSelectionGibbsCI' function. This would involve
// implementing a Markov Chain Monte Carlo (MCMC) algorithm for Bayesian model
// selection, which is a complex statistical task.
struct ModelSelectionGibbsResult {
    arma::mat postSample;
    arma::vec margpp;
    arma::vec postMode;
    double postModeProb;
    arma::vec postProb;
};

ModelSelectionGibbsResult modelSelectionGibbsCI_cpp(
    const arma::vec &postMode_in, double postModeProb_in, int knownphi, int familyint,
    const arma::vec &prior, const arma::vec &priorgr, int niter, int thinning,
    int burnin, int ndeltaini, const arma::ivec &deltaini, const arma::ivec &includevars,
    int n, int p, const arma::vec &ystd, const arma::vec &uncens, double sumy2,
    double sumy, double sumlogyfact, const arma::vec &xstd_vec, const arma::vec &colsumsx,
    bool hasXtX, const arma::mat &XtX, const arma::vec &ytX, const std::string &method,
    const std::string &adj_overdisp, const std::string &hesstype,
    const std::string &optimMethod, int optim_maxit, const arma::vec &thinit,
    bool usethinit, int B, double alpha, double lambda, const arma::vec &phi,
    double tau, double taugroup, double taualpha, bool fixatanhalpha, double r,
    const arma::vec &prDelta, const arma::vec &prDeltap, const arma::vec &parprDeltap,
    const arma::vec &prConstr, const arma::vec &prConstrp,
    const arma::vec &parprConstrp, const arma::vec &groups, int ngroups,
    const arma::vec &nvaringroup, const arma::vec &constraints,
    const arma::vec &invconstraints, int verbose) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'modelSelectionGibbsCI' function's behavior and the MCMC algorithm.
    ModelSelectionGibbsResult result;
    result.postSample.zeros(0, (familyint != 0) ? p : (p + 2));
    result.margpp.zeros((familyint != 0) ? p : (p + 4));
    result.postMode     = postMode_in;
    result.postModeProb = postModeProb_in;
    result.postProb.zeros(1);
    return result;
}

// Placeholder for R's 'listmodels' function.
arma::mat listmodels_cpp(const arma::vec &vars2list, const arma::vec &includevars,
                         const std::vector<arma::vec> &constraints_list,
                         const arma::vec &nvaringroup, int maxvars) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'listmodels' function's behavior for generating model combinations.
    arma::mat models(0, includevars.n_elem);
    return models;
}

// Placeholder for R's 'modelSelectionEnumCI' function. This would involve
// evaluating all possible model combinations, which can be computationally
// expensive.
struct ModelSelectionEnumResult {
    arma::vec postMode;
    double postModeProb;
    arma::vec postProb;
};

ModelSelectionEnumResult modelSelectionEnumCI_cpp(
    int nmodels, const arma::imat &models_in, int knownphi, int familyint,
    const arma::vec &prior, const arma::vec &priorgr, int n, int p,
    const arma::vec &ystd, const arma::vec &uncens, double sumy2, double sumy,
    double sumlogyfact, const arma::vec &xstd_vec, const arma::vec &colsumsx,
    bool hasXtX, const arma::mat &XtX, const arma::vec &ytX, const std::string &method,
    const std::string &adj_overdisp, const std::string &hesstype,
    const std::string &optimMethod, int optim_maxit, const arma::vec &thinit,
    bool usethinit, int B, double alpha, double lambda, const arma::vec &phi,
    double tau, double taugroup, double taualpha, bool fixatanhalpha, double r,
    const arma::vec &prDelta, const arma::vec &prDeltap, const arma::vec &parprDeltap,
    const arma::vec &prConstr, const arma::vec &prConstrp,
    const arma::vec &parprConstrp, const arma::vec &groups, int ngroups,
    const arma::vec &nvaringroup, const arma::vec &constraints,
    const arma::vec &invconstraints, int verbose) {
    // This is a placeholder. The actual implementation would depend on the
    // R 'modelSelectionEnumCI' function's behavior for evaluating models.
    ModelSelectionEnumResult result;
    result.postMode.zeros((familyint != 0) ? p : (p + 2));
    result.postModeProb = 0.0;
    result.postProb.zeros(nmodels);
    return result;
}

// Main C++ function equivalent to the R function
struct MSFitResult {
    arma::mat postSample;
    arma::vec margpp;
    arma::vec postMode;
    double postModeProb;
    arma::vec postProb;
    std::vector<std::string> modelid;
    arma::vec postmean;
    arma::vec postvar;
    std::string family;
    int p;
    bool enumerate;
    std::map<std::string, arma::vec> priors;
    arma::vec ystd;
    arma::mat xstd;
    arma::vec groups;
    arma::vec constraints;
    arma::mat stdconstants;
    std::string outcometype;
    std::map<std::string, std::string> call;
    arma::mat models; // Only if enumerate is true
};

MSFitResult msfit_cpp(
    const arma::vec &y, const arma::mat &x, const arma::mat &data,
    const arma::mat &smoothterms, int nknots, const arma::vec &groups_in,
    const arma::vec &constraints_in, bool center, bool scale, bool enumerate_in,
    const arma::vec &includevars_in, const arma::mat &models_in, int maxvars,
    int niter, int thinning, int burnin, const std::string &family,
    const arma::vec &priorCoef_in, const arma::vec &priorGroup_in,
    const arma::vec &priorDelta_in, const arma::vec &priorConstraints_in,
    const arma::vec &priorVar_in, const arma::vec &priorSkew_in,
    const arma::vec &phi_in, const arma::vec &deltaini_in,
    const std::string &initSearch, const std::string &method,
    const std::string &adj_overdisp, const std::string &hess,
    const std::string &optimMethod, int optim_maxit, const std::string &initpar,
    int B, bool XtXprecomp, bool verbose) {
    // tmp <- formatInputdata(y = y, x = x, data = data, smoothterms = smoothterms,
    //    nknots = nknots, family = family)
    FormattedData tmp      = formatInputdata_cpp(y, x, data, smoothterms, nknots, family);
    arma::mat xstd         = tmp.x;
    arma::vec ystd         = tmp.y;
    std::string formula    = tmp.formula;
    arma::vec splineDegree = tmp.splineDegree;

    // if (!is.null(tmp$groups))
    //    groups <- tmp$groups
    arma::vec groups = (tmp.groups.n_elem > 0) ? tmp.groups : groups_in;

    // if (length(groups) != ncol(x))
    //    stop(paste("groups has the wrong length. It should have length",
    //        ncol(x)))
    if (groups.n_elem != x.n_cols) {
        std::cerr << "Error: groups has the wrong length. It should have length "
                << x.n_cols << std::endl;
        exit(1);
    }
    bool hasgroups     = tmp.hasgroups;
    arma::vec isgroups = tmp.isgroups;

    // if (!is.null(tmp$constraints))
    //    constraints <- tmp$constraints
    arma::vec constraints = (tmp.constraints.n_elem > 0) ? tmp.constraints : constraints_in;

    std::string outcometype            = tmp.outcometype;
    arma::vec uncens                   = tmp.uncens;
    arma::vec ordery                   = tmp.ordery;
    std::vector<std::string> typeofvar = tmp.typeofvar;

    // call <- list(formula = formula, smoothterms = NULL, splineDegree = splineDegree,
    //    nknots = nknots)
    std::map<std::string, std::string> call_cpp;
    call_cpp["formula"]      = formula;
    call_cpp["smoothterms"]  = "NULL"; // In C++, we can just not include it if missing
    call_cpp["splineDegree"] = arma::conv_to<std::string>::from(splineDegree); // Needs proper conversion
    call_cpp["nknots"]       = std::to_string(nknots);

    // if (!missing(smoothterms))
    //    call$smoothterms <- smoothterms
    // In C++, we would check if 'smoothterms' was provided as an argument.

    int p = x.n_cols;
    int n = y.n_elem;

    // if (is.numeric(includevars)) {
    //    tmp = rep(FALSE, p)
    //    if (max(includevars) > p)
    //        stop(paste("includevars contains index ", max(includevars),
    //            " but the design matrix only has ", p, " columns",
    //            sep = ""))
    //    tmp[includevars] = TRUE
    //    includevars = tmp
    // }
    arma::ivec includevars_cpp = arma::zeros<arma::ivec>(p);
    if (includevars_in.n_elem > 0) {
        if (arma::max(includevars_in) > p) {
            std::cerr << "Error: includevars contains index " << arma::max(includevars_in)
                    << " but the design matrix only has " << p << " columns" << std::endl;
            exit(1);
        }
        for (arma::uword i = 0; i < includevars_in.n_elem; ++i) {
            if (includevars_in(i) > 0 && includevars_in(i) <= p) {
                includevars_cpp(includevars_in(i) - 1) = 1;
            }
        }
    }

    // if (length(includevars) != ncol(x) | (!is.logical(includevars)))
    //    stop("includevars must be a logical vector of length ncol(x)")
    if (includevars_cpp.n_elem != x.n_cols) {
        std::cerr << "Error: includevars must be a logical vector of length ncol(x)" << std::endl;
        exit(1);
    }

    // if (missing(maxvars))
    //    maxvars = ifelse(family == "auto", p + 2, p)
    if (maxvars == -1) {
        // Assuming -1 or some special value indicates missing
        maxvars = (family == "auto") ? (p + 2) : p;
    }

    // if (maxvars <= sum(includevars))
    //    stop("maxvars must be >= sum(includevars)")
    if (maxvars <= arma::sum(includevars_cpp)) {
        std::cerr << "Error: maxvars must be >= sum(includevars)" << std::endl;
        exit(1);
    }

    // if (missing(priorCoef)) {
    //    defaultprior = defaultmom(outcometype = outcometype,
    //        family = family)
    //    priorCoef = defaultprior$priorCoef
    //    priorVar = defaultprior$priorVar
    // }
    arma::vec priorCoef = priorCoef_in;
    arma::vec priorVar  = priorVar_in;
    if (priorCoef.is_empty()) {
        DefaultPrior defaultprior = defaultmom_cpp(outcometype, family);
        priorCoef                 = defaultprior.priorCoef;
        priorVar                  = defaultprior.priorVar;
    }

    // if (missing(priorGroup)) {
    //    if (length(groups) == length(unique(groups))) {
    //        priorGroup = priorCoef
    //    }
    //    else {
    //        priorGroup = groupzellnerprior(tau = n)
    //    }
    // }
    arma::vec priorGroup = priorGroup_in;
    if (priorGroup.is_empty()) {
        arma::vec unique_groups = arma::unique(groups);
        if (groups.n_elem == unique_groups.n_elem) {
            priorGroup = priorCoef;
        } else {
            priorGroup = groupzellnerprior_cpp(n);
        }
    }

    // tmp = codeGroupsAndConstraints(p = p, groups = groups, constraints = constraints)
    GroupConstraintInfo group_info = codeGroupsAndConstraints_cpp(p, groups, constraints);
    int ngroups                    = group_info.ngroups;
    constraints                    = group_info.constraints;
    arma::vec invconstraints       = group_info.invconstraints;
    arma::vec nvaringroup          = group_info.nvaringroup;
    groups                         = group_info.groups;

    // if (missing(models)) {
    //    if (missing(enumerate))
    //        enumerate = ifelse(ngroups < 15, TRUE, FALSE)
    // }
    bool enumerate_cpp   = enumerate_in;
    arma::mat models_cpp = models_in;
    if (models_cpp.n_elem == 0) {
        if (enumerate_in == -1) {
            // Assuming -1 indicates missing
            enumerate_cpp = (ngroups < 15);
        }
    } else {
        enumerate_cpp = true;
    }

    // else {
    //    enumerate = TRUE
    // }

    // if (!is.vector(y)) {
    //    y <- as.double(as.vector(y))
    // }
    // else {
    //    y <- as.double(y)
    // }
    // Handled by using arma::vec for y

    // if (!is.matrix(x))
    //    x <- as.matrix(x)
    // Handled by using arma::mat for x

    // mx = colMeans(x)
    arma::rowvec mx = arma::mean(x, 0);
    // sx = sqrt(colMeans(x^2) - mx^2) * sqrt(n/(n - 1))
    arma::rowvec sx = arma::sqrt(arma::mean(arma::pow(x, 2), 0) - arma::pow(mx, 2));
    if (n > 1) {
        sx = sx * std::sqrt(static_cast<double>(n) / (n - 1));
    } else {
        sx.fill(1.0); // Handle case n=1 to avoid division by zero
    }
    // ct = (sx == 0)
    arma::urowvec ct = (sx == 0);

    // if (any(is.na(ct)))
    //    stop("x contains NAs, this is currently not supported, please remove the NAs")
    if (arma::any(arma::is_nan(x))) {
        std::cerr << "Error: x contains NAs, this is currently not supported, please remove the NAs" << std::endl;
        exit(1);
    }

    // if (sum(ct) > 1)
    //    stop("There are >1 constant columns in x (e.g. two intercepts)")
    if (arma::sum(arma::conv_to<arma::vec>::from(ct)) > 1) {
        std::cerr << "Error: There are >1 constant columns in x (e.g. two intercepts)" << std::endl;
        exit(1);
    }

    double my                = 0.0;
    arma::rowvec mx_centered = arma::zeros<arma::rowvec>(p);
    // if (!center) {
    //    my = 0
    //    mx = rep(0, p)
    // }
    if (center) {
        // else {
        //    my = mean(y)
        my          = arma::mean(y);
        mx_centered = mx;
        // }
    }

    double sy              = 1.0;
    arma::rowvec sx_scaled = arma::ones<arma::rowvec>(p);
    // if (!scale) {
    //    sy = 1
    //    sx = rep(1, p)
    // }
    if (scale) {
        // else {
        //    sy = sd(y)
        sy = arma::stddev(y, 0);
        // }
        sx_scaled = sx;
    }

    // mx[typeofvar == "factor"] = 0
    // sx[typeofvar == "factor"] = 1
    for (arma::uword i = 0; i < typeofvar.size(); ++i) {
        if (typeofvar[i] == "factor") {
            mx_centered(i) = 0.0;
            sx_scaled(i)   = 1.0;
        }
    }

    // if (!(outcometype %in% c("Continuous", "Survival"))) {
    //    my = 0
    //    sy = 1
    // }
    if (outcometype != "Continuous" && outcometype != "Survival") {
        my = 0.0;
        sy = 1.0;
    }

    // ystd = (y - my)/sy
    arma::vec ystd_centered_scaled = (y - my) / sy;
    // xstd = x
    // xstd[, !ct] = t((t(x[, !ct]) - mx[!ct])/sx[!ct])
    arma::mat xstd_centered_scaled = x;
    arma::uvec not_ct_indices      = arma::find(ct == 0);
    for (arma::uword i = 0; i < not_ct_indices.n_elem; ++i) {
        arma::uword col_index               = not_ct_indices(i);
        xstd_centered_scaled.col(col_index) = (xstd_centered_scaled.col(col_index) - mx_centered(col_index)) /
                                              sx_scaled(col_index);
    }

    int knownphi = 0;
    arma::vec phi;
    // if (missing(phi)) {
    //    knownphi <- as.integer(0)
    //    phi <- double(0)
    // }
    if (phi_in.n_elem > 0) {
        // else {
        //    knownphi <- as.integer(1)
        knownphi = 1;
        //    phi <- as.double(phi)
        phi = phi_in;
        // }
    } else {
        phi.zeros(0);
    }

    // stdconstants = rbind(c(my, sy), cbind(mx, sx))
    arma::mat stdconstants(2, 2 + p, arma::fill::zeros);
    stdconstants(0, 0) = my;
    stdconstants(0, 1) = sy;
    stdconstants.row(1).subvec(0, 1).zeros(); // Ensure these are 0 as per R code structure
    stdconstants.row(1).subvec(2, 2 + p - 1) = mx;
    stdconstants.row(0).subvec(2, 2 + p - 1) = sx;
    // colnames(stdconstants) = c("shift", "scale") - Not directly applicable in C++

    arma::ivec deltaini_cpp;
    int ndeltaini;
    // if (missing(deltaini)) {
    //    deltaini = as.integer(which(includevars) - 1)
    //    ndeltaini = as.integer(length(deltaini))
    // }
    if (deltaini_in.n_elem == 0) {
        arma::uvec include_indices = arma::find(includevars_cpp == 1);
        deltaini_cpp.set_size(include_indices.n_elem);
        for (arma::uword i = 0; i < include_indices.n_elem; ++i) {
            deltaini_cpp(i) = include_indices(i);
        }
        ndeltaini = deltaini_cpp.n_elem;
    } else {
        // else {
        //    if (length(deltaini) != p)
        //        stop("deltaini must be of length ncol(x)")
        if (deltaini_in.n_elem != p) {
            std::cerr << "Error: deltaini must be of length ncol(x)" << std::endl;
            exit(1);
        }
        //    if (!is.logical(deltaini)) {
        //        stop("deltaini must be of type logical")
        //    }
        // Assuming deltaini_in is a logical vector represented by 0s and 1s
        if (arma::any(deltaini_in != 0 && deltaini_in != 1)) {
            std::cerr << "Error: deltaini must be of type logical (0 or 1)" << std::endl;
            exit(1);
        }
        //    else {
        //        ndeltaini <- as.integer(sum(deltaini | includevars))
        arma::ivec combined_delta = arma::conv_to<arma::ivec>::from(
            arma::max(arma::join_cols(deltaini_in, arma::conv_to<arma::vec>::from(includevars_cpp)), 1));
        ndeltaini = arma::sum(combined_delta);
        //        deltaini <- as.integer(which(deltaini | includevars) -
        arma::uvec combined_indices = arma::find(combined_delta == 1);
        deltaini_cpp.set_size(combined_indices.n_elem);
        for (arma::uword i = 0; i < combined_indices.n_elem; ++i) {
            deltaini_cpp(i) = combined_indices(i);
        }
        //            1)
        //    }
        // }
    }

    // thinit = getthinit(y = y, x = xstd, family = family, initpar = initpar,
    //    enumerate = enumerate)
    ThinitInfo thinit_info = getthinit_cpp(y, xstd_centered_scaled, family, initpar, enumerate_cpp);
    bool usethinit         = thinit_info.usethinit;
    arma::vec thinit       = thinit_info.thinit;

    // method <- formatmsMethod(method = method, usethinit = usethinit,
    //    optimMethod = optimMethod, optim_maxit = optim_maxit,
    //    priorCoef = priorCoef, priorGroup = priorGroup, knownphi = knownphi,
    //    outcometype = outcometype, family = family, hasgroups = hasgroups,
    //    adj.overdisp = adj.overdisp, hess = hess)
    MethodInfo method_info = formatmsMethod_cpp(method, usethinit, optimMethod, optim_maxit,
                                                priorCoef, priorGroup, knownphi, outcometype,
                                                family, hasgroups, adj_overdisp, hess);
    std::string optimMethod_cpp  = method_info.optimMethod;
    int optim_maxit_cpp          = method_info.optim_maxit;
    std::string adj_overdisp_cpp = method_info.adj_overdisp;
    std::string hesstype_cpp     = method_info.hesstype;
    std::string method_cpp       = method_info.method;

    int niter_cpp    = niter;
    int burnin_cpp   = burnin;
    int thinning_cpp = thinning;
    int B_cpp        = B;

    double sumy2          = arma::sum(arma::pow(ystd_centered_scaled, 2));
    double sumy           = arma::sum(ystd_centered_scaled);
    arma::rowvec ytX      = ystd_centered_scaled.t() * xstd_centered_scaled;
    arma::rowvec colsumsx = arma::sum(xstd_centered_scaled, 0);

    arma::mat XtX;
    bool hasXtX_cpp;
    // if (XtXprecomp) {
    //    XtX = t(xstd) %*% xstd
    //    hasXtX = as.logical(TRUE)
    // }
    if (XtXprecomp) {
        XtX        = xstd_centered_scaled.t() * xstd_centered_scaled;
        hasXtX_cpp = true;
    } else {
        // else {
        //    XtX = double(0)
        XtX.zeros(0, 0);
        //    hasXtX = as.logical(FALSE)
        hasXtX_cpp = false;
        // }
    }

    // ffamily = formatFamily(family, issurvival = length(uncens) >
    //    0)
    FamilyInfo family_info = formatFamily_cpp(family, uncens.n_elem > 0);
    int familyint          = family_info.familyint;
    int familygreedy       = family_info.familygreedy;

    double sumlogyfact = 0.0;
    // if (familyint == 22) {
    //    sumlogyfact = as.double(sum(lgamma(ystd + 1))) // lgamma not in armadillo
    // }
    // else {
    //    sumlogyfact = as.double(0)
    // }
    // Assuming familyint == 22 corresponds to Poisson or similar where this might be relevant.
    // This part requires careful translation based on the actual R implementation.

    arma::vec nn_vec;
    // if (!is.null(colnames(xstd))) {
    //    nn <- colnames(xstd) // Column names not directly available in Armadillo
    // }
    // else {
    //    nn <- paste("x", 1:ncol(xstd), sep = "")
    nn_vec.set_size(xstd_centered_scaled.n_cols);
    for (arma::uword i = 0; i < xstd_centered_scaled.n_cols; ++i) {
        nn_vec(i) = i + 1; // Placeholder for column names
    }

    // tmp = formatmsPriorsMarg(priorCoef = priorCoef, priorGroup = priorGroup,
    //    priorVar = priorVar, priorSkew = priorSkew, n = n)
    MarginalPriorInfo marg_prior_info = formatmsPriorsMarg_cpp(priorCoef, priorGroup, priorVar, priorSkew_in, n);
    double r                          = marg_prior_info.r;
    arma::vec prior                   = marg_prior_info.prior;
    arma::vec priorgr                 = marg_prior_info.priorgr;
    double tau                        = marg_prior_info.tau;
    double taugroup                   = marg_prior_info.taugroup;
    double alpha                      = marg_prior_info.alpha;
    double lambda                     = marg_prior_info.lambda;
    double taualpha                   = marg_prior_info.taualpha;
    bool fixatanhalpha                = marg_prior_info.fixatanhalpha;
    priorCoef                         = marg_prior_info.priorCoef;
    priorGroup                        = marg_prior_info.priorGroup;

    // priorConstraints <- defaultpriorConstraints(priorDelta,
    //    priorConstraints)
    arma::vec priorConstraints_cpp = defaultpriorConstraints_cpp(priorDelta_in, priorConstraints_in);

    // tmp = formatmsPriorsModel(priorDelta = priorDelta, priorConstraints = priorConstraints,
    //    constraints = constraints)
    ModelPriorInfo model_prior_info = formatmsPriorsModel_cpp(priorDelta_in, priorConstraints_cpp, constraints);
    arma::vec prDelta               = model_prior_info.prDelta;
    arma::vec prDeltap              = model_prior_info.prDeltap;
    arma::vec parprDeltap           = model_prior_info.parprDeltap;
    arma::vec prConstr              = model_prior_info.prConstr;
    arma::vec prConstrp             = model_prior_info.prConstrp;
    arma::vec parprConstrp          = model_prior_info.parprConstrp;

    arma::vec postMode;
    double postModeProb;
    arma::mat postSample;
    arma::vec margpp;
    arma::vec postProb;
    std::vector<std::string> modelid;
    arma::mat models_out;

    if (!enumerate_cpp) {
        // includevars <- as.integer(includevars)
        arma::ivec includevars_int = arma::conv_to<arma::ivec>::from(includevars_cpp);
        // if (familyint == 0) {
        //    postMode <- rep(as.integer(0), p + 2)
        // }
        // else {
        //    postMode <- rep(as.integer(0), p)
        // }
        postMode.zeros((familyint == 0) ? (p + 2) : p);
        postModeProb = 0.0;
        // postModeProb <- double(1) - Not directly translatable

        // if (initSearch == "greedy") {
        if (initSearch == "greedy") {
            // niterGreed <- as.integer(100)
            int niterGreed = 100;
            // ans = greedyVarSelCI(knownphi, familygreedy,
            //    prior, priorgr, niterGreed, ndeltaini, deltaini, includevars,
            //    n, p, ystd, uncens, sumy2, sumy, sumlogyfact,
            //    xstd, colsumsx, hasXtX, XtX, ytX, method, adj.overdisp,
            //    hesstype, optimMethod, optim_maxit, thinit,
            //    usethinit, B, alpha, lambda, phi, tau, taugroup,
            //    taualpha, fixatanhalpha, r, prDelta, prDeltap,
            //    parprDeltap, prConstr, prConstrp, parprConstrp,
            //    groups, ngroups, nvaringroup, constraints, invconstraints,
            //    as.integer(verbose))
            GreedyVarSelResult greedy_result = greedyVarSelCI_cpp(
                knownphi, familygreedy, prior, priorgr, niterGreed, ndeltaini, deltaini_cpp,
                includevars_int, n, p, ystd_centered_scaled, uncens, sumy2, sumy, sumlogyfact,
                xstd_centered_scaled, arma::conv_to<arma::vec>::from(colsumsx), hasXtX_cpp, XtX,
                arma::conv_to<arma::vec>::from(ytX), method_cpp, adj_overdisp_cpp, hesstype_cpp,
                optimMethod_cpp, optim_maxit_cpp, thinit, usethinit, B_cpp, alpha, lambda,
                phi, tau, taugroup, taualpha, fixatanhalpha, r, prDelta, prDeltap, parprDeltap,
                prConstr, prConstrp, parprConstrp, groups, ngroups, nvaringroup, constraints,
                invconstraints, verbose);
            // postMode <- ans[[1]]
            postMode = greedy_result.postMode;
            // postModeProb <- ans[[2]]
            postModeProb = greedy_result.postModeProb;
            // if (familyint == 0) {
            //    postMode <- as.integer(c(postMode, 0, 0))
            //    postModeProb <- as.double(postModeProb - 2 *
            //        log(2))
            // }
            if (familyint == 0) {
                arma::vec temp_postMode = arma::zeros<arma::vec>(p + 2);
                for (arma::uword i = 0; i < postMode.n_elem; ++i) {
                    temp_postMode(i) = postMode(i);
                }
                postMode = temp_postMode;
                postModeProb -= 2 * std::log(2.0);
            }
            // postMode[includevars == 1] <- TRUE
            for (arma::uword i = 0; i < includevars_cpp.n_elem; ++i) {
                if (includevars_cpp(i) == 1) {
                    postMode(i) = 1;
                }
            }
            // ndeltaini <- as.integer(sum(postMode))
            ndeltaini = arma::sum(postMode);
            // deltaini <- as.integer(which(as.logical(postMode)) -
            arma::uvec postMode_indices = arma::find(postMode == 1);
            deltaini_cpp.set_size(postMode_indices.n_elem);
            for (arma::uword i = 0; i < postMode_indices.n_elem; ++i) {
                deltaini_cpp(i) = postMode_indices(i);
            }
            //    1)
        } else if (initSearch == "SCAD") {
            // else if (initSearch == "SCAD") {
            //    if (verbose)
            //        cat("Initializing via SCAD cross-validation...")
            if (verbose) {
                std::cout << "Initializing via SCAD cross-validation..." << std::endl;
            }
            // deltaini <- rep(TRUE, ncol(xstd))
            arma::ivec deltaini_scad       = arma::ones<arma::ivec>(xstd_centered_scaled.n_cols);
            arma::uvec not_ct_indices_bool = arma::find(arma::conv_to<arma::ivec>::from(!ct) == 1);
            arma::mat xstd_no_ct           = xstd_centered_scaled.cols(not_ct_indices_bool);
            arma::vec ystd_centered        = ystd_centered_scaled - arma::mean(ystd_centered_scaled);
            arma::ivec includevars_no_ct   = arma::conv_to<arma::ivec>::from(includevars_cpp.elem(not_ct_indices_bool));

            arma::ivec scad_init_result = initialize_scad_cpp(xstd_no_ct, ystd_centered, xstd_no_ct.n_cols,
                                                              includevars_no_ct, verbose);
            arma::uvec scad_indices = not_ct_indices_bool(arma::find(scad_init_result == 1));
            deltaini_scad.zeros();
            for (arma::uword i = 0; i < scad_indices.n_elem; ++i) {
                deltaini_scad(scad_indices(i)) = 1;
            }

            // deltaini[!ct] <- ncvreg(...) // Placeholder
            // deltaini[includevars == 1] <- TRUE
            for (arma::uword i = 0; i < includevars_cpp.n_elem; ++i) {
                if (includevars_cpp(i) == 1) {
                    deltaini_scad(i) = 1;
                }
            }
            // ndeltaini <- as.integer(sum(deltaini))
            ndeltaini = arma::sum(deltaini_scad);
            // deltaini <- as.integer(which(deltaini) - 1)
            arma::uvec deltaini_scad_indices = arma::find(deltaini_scad == 1);
            deltaini_cpp.set_size(deltaini_scad_indices.n_elem);
            for (arma::uword i = 0; i < deltaini_scad_indices.n_elem; ++i) {
                deltaini_cpp(i) = deltaini_scad_indices(i);
            }
            // if (verbose)
            //    cat(" Done\n")
        }
        // ans <- modelSelectionGibbsCI(postMode, postModeProb,
        //    knownphi, familyint, prior, priorgr, niter, thinning,
        //    burnin, ndeltaini, deltaini, includevars, n, p,
        //    ystd, uncens, sumy2, sumy, sumlogyfact, as.double(xstd),
        //    colsumsx, hasXtX, XtX, ytX, method, adj.overdisp,
        //    hesstype, optimMethod, optim_maxit, thinit, usethinit,
        //    B, alpha, lambda, phi, tau, taugroup, taualpha,
        //    fixatanhalpha, r, prDelta, prDeltap, parprDeltap,
        //    prConstr, prConstrp, parprConstrp, groups, ngroups,
        //    nvaringroup, constraints, invconstraints, as.integer(verbose))
        ModelSelectionGibbsResult gibbs_result = modelSelectionGibbsCI_cpp(
            postMode, postModeProb, knownphi, familyint, prior, priorgr, niter_cpp,
            thinning_cpp, burnin_cpp, ndeltaini, deltaini_cpp, includevars_int, n, p,
            ystd_centered_scaled, uncens, sumy2, sumy, sumlogyfact,
            arma::conv_to<arma::vec>::from(xstd_centered_scaled.t()),
            arma::conv_to<arma::vec>::from(colsumsx), hasXtX_cpp, XtX, arma::conv_to<arma::vec>::from(ytX),
            method_cpp, adj_overdisp_cpp, hesstype_cpp, optimMethod_cpp, optim_maxit_cpp,
            thinit, usethinit, B_cpp, alpha, lambda, phi, tau, taugroup, taualpha,
            fixatanhalpha, r, prDelta, prDeltap, parprDeltap, prConstr, prConstrp,
            parprConstrp, groups, ngroups, nvaringroup, constraints, invconstraints, verbose);
        // postSample <- matrix(ans[[1]], ncol = ifelse(familyint !=
        //    0, p, p + 2))
        postSample = gibbs_result.postSample;
        // margpp <- ans[[2]]
        margpp = gibbs_result.margpp;
        // postMode <- ans[[3]]
        postMode = gibbs_result.postMode;
        // postModeProb <- ans[[4]]
        postModeProb = gibbs_result.postModeProb;
        // postProb <- ans[[5]]
        postProb = gibbs_result.postProb;
        // margpp[includevars == 1] = 1
        for (arma::uword i = 0; i < includevars_cpp.n_elem; ++i) {
            if (includevars_cpp(i) == 1) {
                margpp(i) = 1.0;
            }
        }
        // postmean = postvar = NULL
        arma::vec postmean_out;
        arma::vec postvar_out;
        postmean = postmean_out;
        postvar  = postvar_out;
        // modelid = apply(postSample[, 1:ncol(xstd), drop = FALSE] ==
        //    1, 1, function(z) paste(which(z), collapse = ","))
        modelid.resize(postSample.n_rows);
        for (arma::uword i = 0; i < postSample.n_rows; ++i) {
            std::string current_model_id = "";
            for (arma::uword j = 0; j < xstd_centered_scaled.n_cols; ++j) {
                if (postSample(i, j) == 1) {
                    if (!current_model_id.empty()) {
                        current_model_id += ",";
                    }
                    current_model_id += std::to_string(j + 1);
                }
            }
            modelid[i] = current_model_id;
        }
    } else {
        // else {
        //    if (verbose)
        //        cat("Enumerating models...\n")
        if (verbose) {
            std::cout << "Enumerating models..." << std::endl;
        }
        // nincludevars = sum(includevars)
        int nincludevars = arma::sum(includevars_cpp);
        // nvars = ifelse(familyint == 0, ncol(xstd) + 2 - nincludevars,
        //    ncol(xstd) - nincludevars)
        int nvars = (familyint == 0)
                        ? (xstd_centered_scaled.n_cols + 2 - nincludevars)
                        : (xstd_centered_scaled.n_cols - nincludevars);

        arma::vec includeenum;
        // if (familyint == 0) {
        //    includeenum = c(includevars[groups + 1], FALSE,
        //        FALSE)
        // }
        if (familyint == 0) {
            includeenum.set_size(groups.n_elem + 2);
            for (arma::uword i = 0; i < groups.n_elem; ++i) {
                if (groups(i) >= 0 && groups(i) < includevars_cpp.n_elem) {
                    includeenum(i) = includevars_cpp(groups(i));
                } else {
                    includeenum(i) = 0; // Handle out-of-bounds group indices
                }
            }
            includeenum(groups.n_elem)     = 0;
            includeenum(groups.n_elem + 1) = 0;
        } else {
            // else {
            //    includeenum = includevars[groups + 1]
            includeenum.set_size(groups.n_elem);
            for (arma::uword i = 0; i < groups.n_elem; ++i) {
                if (groups(i) >= 0 && groups(i) < includevars_cpp.n_elem) {
                    includeenum(i) = includevars_cpp(groups(i));
                } else {
                    includeenum(i) = 0; // Handle out-of-bounds group indices
                }
            }
            // }
        }

        arma::mat models_to_use = models_cpp;
        // if (missing(models)) {
        //    models = listmodels(vars2list = 1:ngroups, includevars = includeenum,
        //        constraints = sapply(constraints, function(z) z +
        //            1), nvaringroup = nvaringroup, maxvars = maxvars)
        // }
        if (models_cpp.n_elem == 0) {
            std::vector<arma::vec> constraints_list;
            for (arma::uword i = 0; i < constraints.n_elem; ++i) {
                constraints_list.push_back(arma::vec(1, arma::fill::value(constraints(i) + 1)));
            }
            arma::vec vars2list(ngroups);
            for (int i = 0; i < ngroups; ++i) vars2list(i) = i + 1;
            models_to_use = listmodels_cpp(vars2list, includeenum, constraints_list, nvaringroup, maxvars);
        } else {
            // else {
            //    if (!is.logical(models))
            //        stop("models must be a logical matrix")
            if (arma::any(models_cpp != 0 && models_cpp != 1)) {
                std::cerr << "Error: models must be a logical matrix (0 or 1)" << std::endl;
                exit(1);
            }
            //    if (ncol(models) != ncol(xstd))
            //        stop(paste("models has", ncol(models), "but x has",
            //            ncol(xstd), "columns"))
            if (models_cpp.n_cols != xstd_centered_scaled.n_cols) {
                std::cerr << "Error: models has " << models_cpp.n_cols << " but x has "
                        << xstd_centered_scaled.n_cols << " columns" << std::endl;
                exit(1);
            }
            // }
        }

        // if (familyint == 0)
        //    models = rbind(cbind(models, FALSE, FALSE), cbind(models,
        //        FALSE, TRUE), cbind(models, TRUE, FALSE), cbind(models,
        //        TRUE, TRUE))
        arma::imat models_int;
        if (familyint == 0) {
            arma::imat models_bool = arma::conv_to<arma::imat>::from(models_to_use);
            arma::imat false_col   = arma::zeros<arma::imat>(models_bool.n_rows, 1);
            arma::imat true_col    = arma::ones<arma::imat>(models_bool.n_rows, 1);
            arma::imat case1       = arma::join_horiz(models_bool, false_col, false_col);
            arma::imat case2       = arma::join_horiz(models_bool, false_col, true_col);
            arma::imat case3       = arma::join_horiz(models_bool, true_col, false_col);
            arma::imat case4       = arma::join_horiz(models_bool, true_col, true_col);
            models_int             = arma::join_vert(case1, case2, case3, case4);
        } else {
            models_int = arma::conv_to<arma::imat>::from(models_to_use);
        }

        int nmodels                     = models_int.n_rows;
        arma::ivec includevars_int_enum = arma::conv_to<arma::ivec>::from(includevars_cpp);

        // ans = modelSelectionEnumCI(nmodels, models, knownphi,
        //    familyint, prior, priorgr, n, p, ystd, uncens, sumy2,
        //    sumy, sumlogyfact, as.double(xstd), colsumsx, hasXtX,
        //    XtX, ytX, method, adj.overdisp, hesstype, optimMethod,
        //    optim_maxit, thinit, usethinit, B, alpha, lambda,
        //    phi, tau, taugroup, taualpha, fixatanhalpha, r,
        //    prDelta, prDeltap, parprDeltap, prConstr, prConstrp,
        //    parprConstrp, groups, ngroups, nvaringroup, constraints,
        //    invconstraints, as.integer(verbose))
        ModelSelectionEnumResult enum_result = modelSelectionEnumCI_cpp(
            nmodels, models_int, knownphi, familyint, prior, priorgr, n, p,
            ystd_centered_scaled, uncens, sumy2, sumy, sumlogyfact,
            arma::conv_to<arma::vec>::from(xstd_centered_scaled.t()),
            arma::conv_to<arma::vec>::from(colsumsx), hasXtX_cpp, XtX, arma::conv_to<arma::vec>::from(ytX),
            method_cpp, adj_overdisp_cpp, hesstype_cpp, optimMethod_cpp, optim_maxit_cpp,
            thinit, usethinit, B_cpp, alpha, lambda, phi, tau, taugroup, taualpha,
            fixatanhalpha, r, prDelta, prDeltap, parprDeltap, prConstr, prConstrp,
            parprConstrp, groups, ngroups, nvaringroup, constraints, invconstraints, verbose);

        // postMode <- ans[[1]]
        postMode = enum_result.postMode;
        // postModeProb <- ans[[2]]
        postModeProb = enum_result.postModeProb;
        // postProb <- ans[[3]]
        postProb = enum_result.postProb;
        // postSample <- matrix(nrow = 0, ncol = ifelse(familyint !=
        //    0, p, p + 2))
        postSample.zeros(0, (familyint != 0) ? p : (p + 2));
        // models <- matrix(models, nrow = nmodels)
        arma::mat models_enum = arma::conv_to<arma::mat>::from(models_int);
        models_out            = models_enum;
        // pp <- exp(postProb - postModeProb)
        arma::vec pp = arma::exp(postProb - postModeProb);
        // pp <- matrix(pp/sum(pp), ncol = 1)
        pp = pp / arma::sum(pp);
        // margpp <- as.vector(t(models) %*% pp)
        margpp = arma::vectorise(models_enum.t() * pp);
        // modelid = apply(models[, 1:ncol(xstd), drop = FALSE] ==
        //    1, 1, function(z) paste(which(z), collapse = ","))
        modelid.resize(models_enum.n_rows);
        for (arma::uword i = 0; i < models_enum.n_rows; ++i) {
            std::string current_model_id = "";
            for (arma::uword j = 0; j < xstd_centered_scaled.n_cols; ++j) {
                if (models_enum(i, j) == 1) {
                    if (!current_model_id.empty()) {
                        current_model_id += ",";
                    }
                    current_model_id += std::to_string(j + 1);
                }
            }
            modelid[i] = current_model_id;
        }

        if (familyint == 0) {
            // modelfam = models[, ncol(xstd) + 1] + 2 * models[,
            //    ncol(xstd) + 2]
            arma::vec modelfam = models_enum.col(xstd_centered_scaled.n_cols) + 2 * models_enum.col(
                                     xstd_centered_scaled.n_cols + 1);
            // margpp = c(margpp[1:ncol(xstd)], sum(pp[modelfam ==
            //    0]), sum(pp[modelfam == 1]), sum(pp[modelfam ==
            //    2]), sum(pp[modelfam == 3]))
            arma::vec margpp_temp = margpp.subvec(0, xstd_centered_scaled.n_cols - 1);
            margpp_temp           = arma::join_cols(margpp_temp, arma::vec({
                                              arma::sum(pp.elem(arma::find(modelfam == 0))),
                                              arma::sum(pp.elem(arma::find(modelfam == 1))),
                                              arma::sum(pp.elem(arma::find(modelfam == 2))),
                                              arma::sum(pp.elem(arma::find(modelfam == 3)))
                                          }));
            margpp = margpp_temp;

            // modeltxt = ifelse(modelfam == 0, "normal", ifelse(modelfam ==
            //    1, "twopiecenormal", ifelse(modelfam == 2, "laplace",
            //    "twopiecelaplace")))
            std::vector<std::string> modeltxt(modelfam.n_elem);
            for (arma::uword i = 0; i < modelfam.n_elem; ++i) {
                if (modelfam(i) == 0) modeltxt[i] = "normal";
                else if (modelfam(i) == 1) modeltxt[i] = "twopiecenormal";
                else if (modelfam(i) == 2) modeltxt[i] = "laplace";
                else if (modelfam(i) == 3) modeltxt[i] = "twopiecelaplace";
            }
            // models = data.frame(modelid = modelid, family = modeltxt,
            //    pp = pp) - Using a struct for output
            // postmean = postvar = NULL
            arma::vec postmean_out;
            arma::vec postvar_out;
            postmean = postmean_out;
            postvar  = postvar_out;
        } else {
            // else {
            //    models = data.frame(modelid = modelid, family = family,
            //        pp = pp) - Using a struct for output
            //    postmean = postvar = NULL
            arma::vec postmean_out;
            arma::vec postvar_out;
            postmean = postmean_out;
            postvar  = postvar_out;
            // }
        }
        // modelid = models$modelid - Already handled
        // models = models[order(models$pp, decreasing = TRUE),
        //    ] - Needs sorting of modelid and potentially models_out based on pp
        arma::uvec sort_indices = arma::sort_index(pp, "descend");
        arma::vec pp_sorted     = pp.elem(sort_indices);
        std::vector<std::string> modelid_sorted(modelid.size());
        arma::mat models_out_sorted = models_out;
        for (arma::uword i = 0; i < sort_indices.n_elem; ++i) {
            modelid_sorted[i] = modelid[sort_indices(i)];
            if (models_out.n_rows > 0) {
                models_out_sorted.row(i) = models_out.row(sort_indices(i));
            }
        }
        modelid    = modelid_sorted;
        models_out = models_out_sorted;
    }

    // if (familyint != 0) {
    //    colnames(postSample) <- names(postMode) <- names(margpp) <- nn
    // }
    // else {
    //    colnames(postSample) <- names(postMode) <- c(nn, "asymmetry",
    //        "laplace")
    //    names(margpp) <- c(nn, "family.normal", "family.tpnormal",
    //        "family.laplace", "family.tplaplace")
    // }
    // Column and name handling not directly applicable in C++ without further structure

    // priors = list(priorCoef = priorCoef, priorGroup = priorGroup,
    //    priorDelta = priorDelta, priorConstraints = priorConstraints,
    //    priorVar = priorVar, priorSkew = priorSkew)
    std::map<std::string, arma::vec> priors_cpp;
    priors_cpp["priorCoef"]        = priorCoef;
    priors_cpp["priorGroup"]       = priorGroup;
    priors_cpp["priorDelta"]       = priorDelta_in;
    priors_cpp["priorConstraints"] = priorConstraints_cpp;
    priors_cpp["priorVar"]         = priorVar;
    priors_cpp["priorSkew"]        = priorSkew_in;

    // if (length(uncens) > 0) {
    //    ystd[ordery] = ystd
    //    uncens[ordery] = uncens
    //    ystd = Surv(time = ystd, event = uncens) // Surv not in Armadillo
    //    xstd[ordery, ] = xstd
    // }
    arma::vec ystd_final = ystd_centered_scaled;
    arma::mat xstd_final = xstd_centered_scaled;
    if (uncens.n_elem > 0 && ordery.n_elem == ystd_centered_scaled.n_elem) {
        ystd_final.zeros(ystd_centered_scaled.n_elem);
        arma::mat xstd_final_temp = arma::zeros<arma::mat>(xstd_centered_scaled.n_rows, xstd_centered_scaled.n_cols);
        for (arma::uword i = 0; i < ordery.n_elem; ++i) {
            if (ordery(i) >= 1 && ordery(i) <= ystd_centered_scaled.n_elem) {
                ystd_final(i) = ystd_centered_scaled(ordery(i) - 1);
            }
        }
        if (ordery.n_elem == xstd_centered_scaled.n_rows) {
            for (arma::uword i = 0; i < ordery.n_elem; ++i) {
                if (ordery(i) >= 1 && ordery(i) <= xstd_centered_scaled.n_rows) {
                    xstd_final_temp.row(i) = xstd_centered_scaled.row(ordery(i) - 1);
                }
            }
            xstd_final = xstd_final_temp;
        }
        // Handling of 'uncens' and 'Surv' would require a specific survival analysis implementation.
    }

    // names(constraints) = paste("group", 0:(length(constraints) -
    //    1)) - Name handling not directly applicable

    MSFitResult result;
    result.postSample   = postSample;
    result.margpp       = margpp;
    result.postMode     = postMode;
    result.postModeProb = postModeProb;
    result.postProb     = postProb;
    result.modelid      = modelid;
    result.postmean     = arma::vec(); // Placeholder
    result.postvar      = arma::vec(); // Placeholder
    result.family       = family;
    result.p            = xstd_centered_scaled.n_cols;
    result.enumerate    = enumerate_cpp;
    result.priors       = priors_cpp;
    result.ystd         = ystd_final;
    result.xstd         = xstd_final;
    result.groups       = groups;
    result.constraints  = constraints;
    result.stdconstants = stdconstants;
    result.outcometype  = outcometype;
    result.call         = call_cpp;
    result.models       = models_out;

    // ans <- list(postSample = postSample, margpp = margpp, postMode = postMode,
    //    postModeProb = postModeProb, postProb = postProb, modelid = modelid,
    //    postmean = postmean, postvar = postvar, family = family,
    //    p = ncol(xstd), enumerate = enumerate, priors = priors,
    //    ystd = ystd, xstd = xstd, groups = groups, constraints = constraints,
    //    stdconstants = stdconstants, outcometype = outcometype,
    //    call = call)
    // if (enumerate) {
    //    ans$models = models
    // }
    // new("msfit", ans) - Class creation not directly translated

    return result;
}

#endif //MOMBF_BRIDGE_H
