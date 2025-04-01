#include <iostream>
#include <armadillo>
#include "controlled-simulation.h"
#include "mombf-bridge.h"

#define DEBUG_PRINTING // Print some additional debug info. Comment to disable.

int main() {
    /* -------------------------------------------------------------------------- */
    /*                               Variable Setup                               */
    /* -------------------------------------------------------------------------- */

    // This section does general setup of variables / settings - this is only needed
    // for development to test the code. Usually these values should be provided from
    // within R

    /* ----------------------- Setup for data generation ------------------------ */
    // TODO I think not all of the descriptions here are correct?

    // Set seed for random number generation
    std::srand(1);

    int n_sim           = 3;     // Number of sim. observations
    int t_sim           = 10;    // Number of sim. time periods
    int nx              = 3;     // Number of regressors
    bool const_val      = false; // Inclusion of a constant
    bool ife            = false; // Inclusion of indiv. fixed effects
    bool tfe            = false; // Inclusion of time fixed effects
    bool iis            = true;  // Inclusion of indicator saturation
    bool sis            = true;  // Inclusion of stepshift saturation
    double p_outl       = 0.0;   // Probability of outlier in a Series
    double p_step       = 0.0;   // Probability of a stepshift in a Series
    int outl_mean       = 0;     // Mean of size of outlier
    int outl_sd         = 0;     // Variance of size of outlier
    int step_mean       = 10;    // Mean of size of stepshift
    double step_sd      = 0.00;  // Variance of size of stepshift
    int error_sd        = 1;     // Variance of the error
    int pos_outl        = 0;     // Position of outlier
    arma::ivec pos_step = {10};  // Position of stepshift

    /* --------------- Setup for the remaining / gibbs parameters --------------- */
    // TODO descriptions are missing

    int i_index         = 0; // Adjusting to 0-based indexing in C++
    int t_index         = 1; // Adjusting to 0-based indexing in C++
    int y_index         = 2; // Adjusting to 0-based indexing in C++
    long long Ndraw     = 5000;
    long long Nburn     = 500;
    std::string b_prior = "g";
    double lambda_b     = 100.0;
    double c0           = 0.001;
    double C0           = 0.001;
    double va           = 1.0;
    double vb           = 1.0;
    double tau          = 1.0;
    bool geweke         = false;

    /* ----------------------------- Generate Data ------------------------------ */

    // Use the above settings to generate random data
    // SimulationOutput sim_result = contr_sim_breaks(
    //     n_sim,
    //     t_sim,
    //     nx,
    //     iis,
    //     sis,
    //     pos_outl,
    //     pos_step,
    //     const_val,
    //     ife,
    //     tfe,
    //     outl_mean,
    //     step_mean,
    //     error_sd
    // );
    //
    // std::cout << sim_result.data << std::endl;
    /* ------------------- GEWEKE Tests not implemented yet! -------------------- */
    // TODO GEWEKE Tests not implemented yet!

    // if(geweke){
    //     library(extraDistr)
    //     lambda_b = lambda_g = 1.234
    //     Ndraw= 100000
    //     c0=3
    //     C0=1
    //     sis=TRUE
    //     b_prior="g"
    //   }


    /* -------------------------------------------------------------------------- */
    /*                              Data Preparation                              */
    /* -------------------------------------------------------------------------- */

    // Prepare the generated data to get it into a form we need

    /* ---------------------------- Data Extraction ----------------------------- */

    arma::mat data = {
        {1, 1, -23.45266676, -2.30936055, -0.55970183, 1.34953562},
        {1, 2, -7.91095366, -1.46139144, 0.21277444, -0.47043324},
        {1, 3, -10.58970802, -0.06862945, -0.12904591, 1.80418798},
        {1, 4, -4.03629912, -0.23848014, 0.31444212, 0.18225964},
        {1, 5, 20.07450224, 0.54917267, 0.18824194, -2.76781820},
        {1, 6, 14.98383100, 1.04541866, 0.01999001, -1.02430597},
        {1, 7, 4.36945658, -0.37946788, 0.32571192, -1.17450487},
        {1, 8, -8.35794726, -0.22696702, -0.28827137, 1.35956658},
        {1, 9, 10.42538523, 0.09199853, -0.72593770, -1.00533029},
        {1, 10, -2.68719828, -0.85342818, -0.29560190, 1.31471504},
        {1, 11, 12.63427833, -0.50768911, -0.93493273, -0.60211238},
        {1, 12, 15.25046563, 1.12089459, -0.43084369, 0.55911127},
        {1, 13, 20.55473321, 0.58630229, -0.68998967, -0.40238799},
        {1, 14, 7.69004829, -1.27681905, 0.06944209, -0.88038903},
        {1, 15, 0.01908629, -0.70259846, 0.75558385, 0.55643298},
        {1, 16, 17.14391279, 1.69326930, 1.02409623, -0.06464430},
        {1, 17, 16.79348380, -0.27937351, -0.68313724, -1.14558556},
        {1, 18, 18.42294041, 1.23860288, 0.89303727, -1.17465964},
        {1, 19, 15.93351524, -0.27851469, -1.51865248, -0.30728423},
        {1, 20, 3.41723235, 0.50750409, 2.13166805, 0.07164625},
        {2, 1, 2.01367442, -0.52964157, -0.84032731, -0.24943370},
        {2, 2, 0.12292821, 1.24425732, 1.89566249, -0.07214188},
        {2, 3, 14.91203433, 0.74228208, 0.57532928, -2.24319079},
        {2, 4, -5.11301574, -0.46578083, 0.89250538, -0.64141230},
        {2, 5, -2.93105027, 0.75886611, 0.61983879, 0.82792110},
        {2, 6, -1.85769585, 0.45025558, 0.09238917, 0.81542852},
        {2, 7, 12.58391461, 0.61893994, -0.07080554, -1.25784309},
        {2, 8, 6.08416188, 0.43634622, -0.35005450, -0.11806201},
        {2, 9, 4.04203852, -0.06634003, -0.35659074, -0.11128336},
        {2, 10, -5.61360271, -0.24980976, -0.11055986, 0.75599906},
        {2, 11, 6.40934135, 1.33751293, -0.26408725, 0.72592710},
        {2, 12, -11.60791505, -1.09736930, -0.61361462, 1.08042091},
        {2, 13, 4.77833559, 1.24434837, 1.19657744, -0.14257870},
        {2, 14, 3.98187096, 0.82541432, -0.90584227, 1.00333061},
        {2, 15, 4.86749764, 1.11334388, 0.66078066, -0.18836024},
        {2, 16, -2.97304347, 0.17946976, 0.54460379, 0.19905417},
        {2, 17, -1.24581570, -1.49029084, -0.48200822, -1.30089775},
        {2, 18, 10.14556265, 0.93577627, -0.38133499, -0.41726513},
        {2, 19, 8.77631832, 0.17082548, -0.53630861, -0.60742406},
        {2, 20, -5.32484498, -0.27814329, -0.17968005, 0.70029078},
        {3, 1, 6.34106622, 1.03389672, -1.18754295, 1.44397345},
        {3, 2, -0.10247840, -0.52907998, 0.71449362, -1.34239109},
        {3, 3, 0.93901800, 0.30634010, -0.78423320, 0.78976343},
        {3, 4, -8.23697372, -1.01782854, 0.22088162, 0.11148258},
        {3, 5, -9.49952944, -0.33937959, 0.82135615, 0.33237310},
        {3, 6, 3.24726127, -0.15963926, -0.69482875, -0.23291455},
        {3, 7, -8.94658045, -2.23058885, -0.09104611, -1.00355918},
        {3, 8, -18.67803117, -0.02682357, 1.79855827, 1.53172391},
        {3, 9, 18.56838167, 2.10492018, 1.05464263, -1.08119123},
        {3, 10, -12.34777296, -0.23481469, 2.47325089, 0.07734140},
        {3, 11, -19.62430162, -1.70293765, 1.14583713, 0.49988338},
        {3, 12, 4.42893308, 1.42905860, -0.12474232, 0.96249562},
        {3, 13, -4.50378125, -0.77566028, -0.58205656, 0.42390948},
        {3, 14, 7.53966294, 0.81322432, 0.35094792, -0.50595807},
        {3, 15, 14.47285868, 1.32722746, -2.44584316, 1.26806703},
        {3, 16, -1.47646208, -0.35883300, -0.44026205, 0.07695551},
        {3, 17, 6.69367344, 1.42208596, 1.30404222, -0.83074969},
        {3, 18, 9.61231686, 0.06315024, -1.24184881, -0.52511487},
        {3, 19, 10.14149005, 0.44613682, -0.63743304, -0.55301099},
        {3, 20, -0.69086139, -2.07200545, -1.89931838, -0.63308381}
    };

    // arma::mat data = sim_result.data;   // Extract the "actual" data from the sim
    arma::vec y = data.col(y_index); // Extract the y vector from the sim data

    // Extract initial X (excluding y, individual index, and time index)
    // Vector of which columns to exclude {0,1,2}
    arma::ivec exclude_cols = {
        y_index, i_index, t_index
    };

    // Our X matrix should have cols(data) - length(excluded_cols) columns
    int num_include_cols = data.n_cols - exclude_cols.n_elem;

    // Pre allocate the Matrix
    arma::mat X_(data.n_rows, num_include_cols);

    // Helper to keep track of the current insertion index into X
    int current_x_col = 0;

    // Loop over all columns of data. If the columns index is one of the ones in
    // exclude_cols, skip, else take the column and insert it into X at current col
    // index. Increment index.
    for (int j = 0; j < data.n_cols; ++j) {
        bool exclude = false;
        for (int k = 0; k < exclude_cols.n_elem; ++k) {
            if (j == exclude_cols(k)) {
                exclude = true;
                break;
            }
        }
        if (!exclude) {
            X_.col(current_x_col) = data.col(j);
            current_x_col++;
        }
    }
    arma::mat X = X_; // Create a copy for further modifications

    // Determine n and t (number of units and time periods)
    arma::vec unique_n_vals = arma::unique(data.col(i_index));
    arma::vec unique_t_vals = arma::unique(data.col(t_index));
    int n                   = unique_n_vals.n_elem;
    int t                   = unique_t_vals.n_elem;
    long N                  = (n) * t;

    /* ----------------------------- Build X Matrix ----------------------------- */
    // Add columns depending on what to include.

    if (const_val) {
        X.insert_cols(X.n_cols, arma::ones(data.n_rows));
    }

    if (ife) {
        // Individual fixed effects
        arma::mat IFE = kronecker_product(arma::eye(n, n), arma::ones(t, 1));
        // No direct equivalent to R's colnames assignment here, but we can keep track of the meaning
        if (const_val) {
            IFE = IFE.cols(1, IFE.n_cols - 1);
        }
        X.insert_cols(X.n_cols, IFE);
    }

    if (tfe) {
        // Time fixed effects
        arma::mat TFE = kronecker_product(arma::ones(n, 1), arma::eye(t, t));
        // No direct equivalent to R's colnames assignment here
        if (const_val) {
            TFE = TFE.cols(1, TFE.n_cols - 1);
        }
        X.insert_cols(X.n_cols, TFE);
    }

    if (tfe && ife && !const_val) {
        std::cerr <<
                "Warning: Both time and unit fixed effects used.\nDropping first indiv. FE to avoid perfect colinearity"
                << std::endl;
        if (!X_.is_empty()) {
            X.shed_col(X_.n_cols); // Assuming the IFE dummies were added after the original X_
        } else if (const_val) {
            // If no original regressors, and IFE was added, it would be at the beginning
            X.shed_col(0);
        } else if (ife) {
            // If only IFE was added, it would be at the beginning
            X.shed_col(0);
        }
    }


    /* ------------------------- Find Matrix dimensions ------------------------- */

    int p = X.n_cols; // Find matrix dimensions

    // Construct z matrix
    arma::mat z_temp(t, t, arma::fill::ones);
    arma::mat z_lower_tri = arma::trimatl(z_temp);
    arma::mat z = z_lower_tri.cols(2, t - 2); // Remove first two and the last column (adjusting indices for 0-based)
    // z = arma::conv_to<arma::mat>::from(arma::round(z));
    // Convert to integer (Armadillo doesn't have direct integer matrix type for mat)


    // Construct Z matrix
    // TODO I think Arma has no integer matrix type? Does that have performance implications?
    arma::mat Z = kronecker_product(arma::eye(n, n), z);
    Z           = arma::conv_to<arma::mat>::from(arma::round(Z)); // Convert to integer


    int r        = Z.n_cols;
    arma::mat XX = arma::trans(X) * X;
    arma::mat ZZ = arma::trans(Z) * Z;

    /* ---------------------------- Prior Parameters ---------------------------- */

    arma::vec b0;
    if (b_prior == "g") {
        b0.zeros(p);
    } else if (b_prior == "f") {
        // b0 = arma::inv(XX) * X.t() * y;
        // Note: For numerical stability, Armadillo provides a more direct way to solve this system:
        b0 = arma::solve(XX, arma::trans(X) * y);
    }

    arma::mat B0     = XX.i() * lambda_b;
    arma::mat B0_inv = XX / lambda_b;


    /* ---------------------------- Starting Values ----------------------------- */

    arma::vec b_i(p, arma::fill::ones);
    arma::vec g_i(r, arma::fill::zeros);
    arma::vec g_incl_i = g_i; // always non zero for rnlp within Gibbs; last value conditional on inclusion
    double s2_i        = 1.0;

    arma::Col<int> w_i(r, arma::fill::zeros);

    // TODO z_cols will be a sequence from 0 to ncol(z) - 1 compared to Rs 1 to ncol(z) because of 0 based indexing
    //  verify that this is correct
    arma::ivec z_cols_temp(z.n_cols);
    for (int i = 0; i < z.n_cols; ++i) {
        z_cols_temp(i) = i;
    }
    arma::ivec z_cols(z.n_cols * n);
    for (int i = 0; i < n; ++i) {
        z_cols.subvec(i * z.n_cols, (i + 1) * z.n_cols - 1) = z_cols_temp;
    }

    // TODO we should just be able to init this with zeros since w_i is always zero
    arma::vec w_1 = arma::conv_to<arma::vec>::from(z_cols) % w_i; // element-wise multiplication

    /* --------------------- Starting Values for HS-Plug-In --------------------- */
    arma::vec lamb_b(p, arma::fill::ones);
    double tau_b = 1.0;
    arma::vec nu_b(p, arma::fill::ones);
    double xi_b = 1.0;

    /* --------------------------- Prep Calculations ---------------------------- */
    double cN        = c0 + N / 2.0;
    arma::mat XX_inv = arma::inv(XX);
    arma::vec Pr_g(n * (t - 2)); // NA in R can be represented by uninitialized values or a specific value if needed


    /* --------------------------------- Store ---------------------------------- */
    long Nstore = 0;
    if (Nburn >= Ndraw) {
        std::cerr << "Error: The number of burn-in exceeds number of draws" << std::endl;
        // Consider throwing an exception or handling this error appropriately in your program
        return 0; // Or exit the function
    } else {
        Nstore = Ndraw - Nburn;
    }

    /* ---------------------------------- Main ---------------------------------- */
    arma::mat b_store(Nstore, p);
    arma::mat g_store(Nstore, r);
    arma::mat s2_store(Nstore, 1);

    arma::Mat<int> w_store(Nstore, r);


    arma::vec g_i_n0 = g_i;
#ifdef DEBUG_PRINTING
    printf("Simulated Data: \n");
    std::cout << data << '\n';

    printf("Extracted y: \n");
    std::cout << y << std::endl;

    printf("Extracted X: \n");
    std::cout << X << std::endl;

    printf("z: \n");
    std::cout << z << '\n';

    printf("Z: \n");
    std::cout << Z << std::endl;

    printf("XX (X'X): \n");
    std::cout << XX << std::endl;

    printf("ZZ (Z'Z): \n");
    std::cout << ZZ << std::endl;

    printf("Prior parameter b0: \n");
    std::cout << b0 << std::endl;

    printf("B0: \n");
    std::cout << B0 << std::endl;

    printf("B0_inv: \n");
    std::cout << B0_inv << std::endl;

    printf("Starting values:\n");
    printf("b_i: \n");
    std::cout << b_i << std::endl;
    printf("g_i: \n");
    std::cout << g_i << std::endl;
    printf("g_incl_i: \n");
    std::cout << g_incl_i << std::endl;
    printf("s2_i: %f\n", s2_i);
    printf("w_i: \n");
    std::cout << w_i << std::endl;
    printf("z_cols: \n");
    std::cout << z_cols << std::endl;
    printf("w_1: \n");
    std::cout << w_1 << std::endl;

    printf("cN: %f\n", cN);
    printf("XX_inv: \n");
    std::cout << XX_inv << std::endl;
    printf("Pr_g (first 10 elements): \n");
    if (Pr_g.n_elem > 10) {
        std::cout << Pr_g.head(10) << std::endl;
    } else {
        std::cout << Pr_g << std::endl;
    }

    printf("Nstore: %ld\n", Nstore);

    printf("b_store dimensions: %zu x %d\n", b_store.n_rows, b_store.n_cols);
    printf("g_store dimensions: %zu x %d\n", g_store.n_rows, g_store.n_cols);
    printf("s2_store dimensions: %zu x %d\n", s2_store.n_rows, s2_store.n_cols);
    printf("w_store dimensions: %zu x %d\n", w_store.n_rows, w_store.n_cols);
#endif


    /* -------------------------------------------------------------------------- */
    /*                               Gibbs Sampler                                */
    /* -------------------------------------------------------------------------- */
    std::cout << "Starting Gibbs Sampler" << std::endl;


    for (int iter = 1 - Nburn; iter <= Nstore; ++iter) {
        /* ------------------------------ Geweke Test ------------------------------- */
        if (geweke) {
            arma::vec e = arma::randn(N) * std::sqrt(s2_i);
            y           = X * b_i + Z * g_i + e;
        }

        /* --------------------------- draw p(s2|a,b,g,y) --------------------------- */
        double CN = C0 + 0.5 * arma::sum(arma::square(y - X * b_i - Z * g_i));
        s2_i      = 1.0 / arma::randg(arma::distr_param(cN, 1.0 / CN)); // 1.0 / CN -> rate = 1/scale

        // s2_i = sqrt(2);

        /* --------------------------- draw p(b|a,g,s2,y) --------------------------- */
        if (b_prior == "hs") {
            // draw p(xi,nu|b,s2,y) ========.====
            arma::vec nu_b(lamb_b.n_elem);
            for (int j = 0; j < lamb_b.n_elem; ++j) {
                nu_b(j) = 1.0 / arma::randg(arma::distr_param(1.0, 1.0 + 1.0 / lamb_b(j)));
                if (nu_b(j) > 1e+8) nu_b(j) = 1e+8;
                if (nu_b(j) < 1e-8) nu_b(j) = 1e-8;
            }
            xi_b = 1.0 / arma::randg(arma::distr_param(1.0, 1.0 + 1.0 / tau_b));
            if (xi_b > 1e+8) xi_b = 1e+8;
            if (xi_b < 1e-8) xi_b = 1e-8;

            // draw p(tau,lamb|xi,nu,b,s2,y).====
            arma::vec lamb_b_new(nu_b.n_elem);
            for (int j = 0; j < nu_b.n_elem; ++j) {
                lamb_b_new(j) = 1.0 / arma::randg(
                                    arma::distr_param(1.0, 1.0 / nu_b(j) + std::pow(b_i(j), 2) / (2.0 * tau_b * s2_i)));
                if (lamb_b_new(j) > 1e+8) lamb_b_new(j) = 1e+8;
                if (lamb_b_new(j) < 1e-8) lamb_b_new(j) = 1e-8;
            }
            lamb_b = lamb_b_new;

            tau_b = 1.0 / arma::randg(arma::distr_param((p + 1.0) / 2.0,
                                                        1.0 / xi_b + arma::sum(arma::square(b_i) / lamb_b) / (
                                                            2.0 * s2_i)));
            if (tau_b > 1e+8) tau_b = 1e+8;
            if (tau_b < 1e-8) tau_b = 1e-8;

            arma::mat A     = XX + (1.0 / tau_b) * arma::diagmat(1.0 / lamb_b);
            arma::mat A_inv = arma::inv(A);
            arma::mat BN    = s2_i * A_inv;
            arma::vec bN    = A_inv * (arma::trans(X) * (y - Z * g_i));

            arma::mat BN_chol;
            bool chol_success = arma::chol(BN_chol, BN);
            if (chol_success) {
                b_i = bN + BN_chol.t() * arma::randn(p);
            } else {
                std::cerr << "Warning: Cholesky decomposition failed. Using multivariate normal from another method." <<
                        std::endl;
                b_i = arma::mvnrnd(bN, BN);
            }
        } else if (b_prior == "g" || b_prior == "f") {
            arma::mat BN_inv = (1.0 / s2_i + 1.0 / lambda_b) * XX;
            arma::mat BN     = (s2_i * lambda_b / (s2_i + lambda_b)) * XX_inv;
            arma::vec bN     = (1.0 / (s2_i + lambda_b)) * (
                               lambda_b * XX_inv * (arma::trans(X) * (y - Z * g_i)) + s2_i * b0);

            arma::mat BN_chol;
            bool chol_success = arma::chol(BN_chol, BN);
            if (chol_success) {
                b_i = bN + BN_chol.t() * arma::randn(p);
            } else {
                std::cerr <<
                        "Warning: Cholesky decomposition failed (g/f prior). Using multivariate normal from another method."
                        << std::endl;
                b_i = arma::mvnrnd(bN, BN);
            }
        }

        // arma::vec fixed_random = {1, -1, 2};
        // b_i                    = fixed_random;

        // ==================== draw p(w|a,b,s2,y) ====================
        arma::vec y_hat = y - X * b_i;

        /* ---------------------------- Model Selection ----------------------------- */

        int nn = 2;

        arma::vec thinit;
        MombfBridge::InitparType initpar_type;
        if (iter == 1 - Nburn) {
            initpar_type = MombfBridge::Auto;
        } else {
            thinit = g_i;
        }

        std::cout << y_hat << std::endl;

        arma::Col<int> post_sample = MombfBridge::modelSelection(y_hat,
                                                                 Z,
                                                                 nn,
                                                                 1,
                                                                 nn - 1,
                                                                 w_i,
                                                                 false,
                                                                 false,
                                                                 true,
                                                                 s2_i,
                                                                 tau,
                                                                 0.348,
                                                                 0.5,
                                                                 thinit,
                                                                 initpar_type
        );


        w_i = post_sample;

        // arma::vec w_i_mod_postSample(r);
        // if (iter == 1 - Nburn) {
        //     w_i_mod_postSample.zeros(); // Placeholder initialization
        // } else {
        //     w_i_mod_postSample.ones(); // Another placeholder
        // }
        // w_i = w_i_mod_postSample;


        if (geweke) {
            if (iter == (1 - Nburn)) {
                std::cout << "\nCAREFUL: No Variable Selection in geweke test!\n" << std::endl;
            }
            w_i.ones(r);
        }


        // ==================== draw p(g|w,a,b,s2,y) ====================
        if (arma::sum(w_i) > 0) {
            arma::vec ans = arma::vec().zeros(Z.n_cols);

            int nonZeroCount = arma::accu(w_i != 0);

            // Create a new matrix with the same number of rows as Z and columns equal to nonZeroCount
            arma::dmat filteredZ(Z.n_rows, nonZeroCount);

            // arma::vec filteredGi(nonZeroCount);

            int currentCol = 0;
            for (int i = 0; i < w_i.size(); ++i) {
                if (w_i(i) != 0) {
                    // Copy the corresponding column from Z to filteredZ
                    filteredZ.col(currentCol) = Z.col(i);
                    // filteredGi(currentCol)    = g_i_n0(i);
                    ++currentCol;
                }
            }

            // Filtered Z -> Nur die spalten wo w_i == 1


            rnlpPost_lm(ans.memptr(),
                        5,
                        4,
                        1,
                        y_hat.memptr(),
                        filteredZ.memptr(),
                        y_hat.size(),
                        filteredZ.n_cols,
                        1,
                        tau_b,
                        c0,
                        C0,
                        1);


            // arma::vec g_i = arma::vec().zeros(w_i.n_elem);
            //
            // currentCol = 0;
            // for (int i = 0; i < w_i.n_elem; ++i) {
            //     g_i[currentCol]      = ans[i];
            //     g_incl_i[currentCol] = ans[i];
            //     currentCol++;
            //     if (w_i(i) == 0) {
            //         currentCol++;
            //     }
            // }

            arma::vec g_i = arma::vec().zeros(w_i.n_elem);

            int ans_pos = 0; // Position tracker for the ans vector
            for (int i = 0; i < w_i.n_elem; ++i) {
                if (w_i(i) == 1) {
                    // When w_i is 1, use the next value from ans
                    g_i(i)      = ans[ans_pos];
                    g_incl_i(i) = ans[ans_pos];
                    ans_pos++; // Only increment ans_pos when we use a value
                }
                // When w_i is 0, g_i stays 0 (already initialized that way)
            }
        } else {
            g_i.zeros(r);
        }

        // TODO
        g_i.zeros(r);
        // g_i(7) = 10;

        // TODO is this a correct nan check?
        if (!arma::find_nan(g_i).is_empty()) {
            std::cerr << "NaNs Produced in g_i at iteration " << iter << std::endl;
        }

        // ==================== store and adjust tau ====================
        if (iter > 0) {
            // R's loop goes from 1-Nburn, so storing starts after burn-in
            b_store.row(iter - 1) = b_i.t();
            g_store.row(iter - 1) = g_i.t();
            w_store.row(iter - 1) = w_i.t();
            s2_store(iter - 1, 0) = s2_i;
        } else {
            // Adjust tau during burn-in (placeholder logic)
            // if (arma::sum(g_i) < n * t * 0.05) {
            //     tau = tau * 5.0 - arma::sum(g_i);
            // }
        }
    }
    std::cout << "Gibbs Sampler finished." << std::endl;

    /* -------------------------------------------------------------------------- */
    /*                        Calculate and Print w_store Means                      */
    /* -------------------------------------------------------------------------- */
    std::cout << "Calculating column means of w_store..." << std::endl;

    // Calculate column means of w_store
    arma::rowvec w_store_means(r);
    for (int j = 0; j < r; ++j) {
        double sum = 0.0;
        for (int i = 0; i < Nstore; ++i) {
            sum += w_store(i, j);
        }
        w_store_means(j) = sum / Nstore;
    }

    // Calculate column means of w_store
    arma::rowvec b_store_means(p);
    for (int j = 0; j < p; ++j) {
        double sum = 0.0;
        for (int i = 0; i < Nstore; ++i) {
            sum += b_store(i, j);
        }
        b_store_means(j) = sum / Nstore;
    }

    arma::rowvec s2_store_means(1);
    double sum = 0.0;
    for (int i = 0; i < Nstore; ++i) {
        sum += s2_store(i, 0);
    }
    s2_store_means(0) = sum / Nstore;


    // Print the results
    std::cout << "Column means of w_store:" << std::endl;
    std::cout << w_store_means.t() << std::endl;


    // Print the results
    std::cout << "Column means of b_store:" << std::endl;
    std::cout << b_store_means.t() << std::endl;

    // Print the results
    std::cout << "Column means of s2_store:" << std::endl;
    std::cout << s2_store_means.t() << std::endl;

    return 0;
}


