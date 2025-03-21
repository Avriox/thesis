#include <iostream>
#include <armadillo>
#include "controlled-simulation.h"

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
    int step_mean       = 5;     // Mean of size of stepshift
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
    SimulationOutput sim_result = contr_sim_breaks(
        n_sim, t_sim, nx, iis, sis, pos_outl, pos_step, const_val, ife, tfe, outl_mean, step_mean, error_sd
    );

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

    arma::mat data = sim_result.data;   // Extract the "actual" data from the sim
    arma::vec y    = data.col(y_index); // Extract the y vector from the sim data

    // Extract initial X (excluding y, individual index, and time index)
    // Vector of which columns to exclude {0,1,2}
    arma::uvec exclude_cols = {
        static_cast<arma::uword>(y_index), static_cast<arma::uword>(i_index), static_cast<arma::uword>(t_index)
    };

    // Our X matrix should have cols(data) - length(excluded_cols) columns

    int num_include_cols = data.n_cols - exclude_cols.n_elem;

    // Pre allocate the Matrix 
    arma::mat X_(data.n_rows, num_include_cols);

    // Helper to keep track of the current insertion index into X
    arma::uword current_x_col = 0;

    // Loop over all columns of data. If the columns index is one of the ones in
    // exclude_cols, skip, else take the column and insert it into X at current col
    // index. Increment index.
    for (arma::uword j = 0; j < data.n_cols; ++j) {
        bool exclude = false;
        for (arma::uword k = 0; k < exclude_cols.n_elem; ++k) {
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

#ifdef DEBUG_PRINTING
    printf("Simulated Data: \n");
    std::cout << data << '\n';

    printf("Extracted y: \n");
    std::cout << y << std::endl;

    printf("Extracted X: \n");
    std::cout << X << std::endl;

#endif

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


#ifdef DEBUG_PRINTING
    printf("z: \n");
    std::cout << z << '\n';

    printf("Z: \n");
    std::cout << Z << std::endl;

#endif


    return 0;
}


