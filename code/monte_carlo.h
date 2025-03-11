//
// Created by Jakob Goldmann on 05.03.25.
//

#ifndef MONTE_CARLO_H
#define MONTE_CARLO_H

#include <random>

// Function to estimate the integral using Monte Carlo simulation in C++ with Armadillo
inline double solve_integral_monte_carlo_cpp(const arma::mat& Ck, const arma::mat& Xk, const arma::vec& y, double sigma_squared, int r, int num_samples) {
    // Function arguments:
    // Ck:  Positive definite symmetric matrix Ck (Armadillo mat object). 'const mat&' means it's passed by reference and won't be modified.
    // Xk:  Matrix Xk (Armadillo mat object).
    // y:   Vector y (Armadillo vec object).
    // sigma_squared: Scalar sigma^2 (double).
    // r:   Exponent r (integer).
    // num_samples: Number of Monte Carlo samples (integer).
    // Returns: Estimated value of the integral (double).

    // ------------------- Step 1: Compute beta_tk and Sigma -------------------

    // 1a. Compute the inverse of Ck.
    arma::mat Gamma = inv(Ck); // arma::inv() function computes the inverse of a matrix.
                          // 'inv(Ck)' calculates Ck^(-1) and stores it in 'Gamma'.

    // 1b. Compute beta_tk = Ck^(-1) * Xk.transpose() * y
    arma::vec beta_tk = Gamma * trans(Xk) * y; // arma::trans(Xk) gets the transpose of Xk (Xk^T).
                                        // '*' operator performs matrix multiplication in Armadillo.
                                        // Order of operations: (Gamma * (trans(Xk) * y)).

    // 1c. Compute the covariance matrix Sigma = sigma^2 * Ck^(-1) = sigma^2 * Gamma
    arma::mat Sigma = sigma_squared * Gamma; // Scalar multiplication with a matrix in Armadillo.


    // ------------------- Step 2: Cholesky decomposition of Sigma -------------------

    arma::mat L = chol(Sigma, "lower"); // arma::chol() function computes the Cholesky decomposition.
                                  // The second argument "lower" specifies to return the lower triangular matrix L.
                                  // Sigma = L * L.t() where L.t() is the transpose of L.


    // ------------------- Step 3: Monte Carlo Simulation -------------------
    // Initialize sum of product terms to 0.0
    double sum_Y = 0.0;

    // Random number generation setup using C++11 <random> library.
    // Create a random number generator engine (e.g., std::mt19937 - Mersenne Twister engine)
    std::random_device rd{}; // Use a hardware entropy source if available, or a pseudo-random source if not.
    std::mt19937 gen{rd()};  // Seed the Mersenne Twister generator with entropy from random_device.
    std::normal_distribution<> d{0, 1}; // Define a normal distribution with mean 0 and standard deviation 1.


    // Loop for the specified number of samples
    for (int i = 0; i < num_samples; ++i) {
        // ------------------- Step 3a: Generate a sample beta from MVN(beta_tk, Sigma) -------------------

        // 3a.i. Generate k independent standard normal random variables (vector z).
        int k = Ck.n_rows; // Get the dimension k from the number of rows of Ck (or Sigma, L, beta_tk size).
        arma::vec z(k);         // Create an Armadillo vector 'z' of size k to store standard normal variables.
        for (int j = 0; j < k; ++j) {
            z(j) = d(gen); // Generate a standard normal sample using the distribution 'd' and generator 'gen', and assign it to z(j).
        }

        // 3a.ii. Compute beta_sample = beta_tk + L * z
        arma::vec beta_sample = beta_tk + L * z; // Matrix-vector multiplication (L * z) and vector addition.


        // ------------------- Step 3b: Calculate the product term Y = prod(beta_j^(2r)) -------------------

        double product_term = 1.0; // Initialize product term for this sample to 1.
        for (int j = 0; j < k; ++j) {
            product_term *= pow(beta_sample(j), 2.0 * r); // Calculate beta_j^(2r) using pow() and multiply to product_term.
                                                        // beta_sample(j) accesses the j-th element of the beta_sample vector (0-indexed).
        }

        // ------------------- Step 3c: Accumulate the product term -------------------
        sum_Y += product_term; // Add the calculated product_term for this sample to the running sum sum_Y.
    }


    // ------------------- Step 4: Calculate the estimated expectation E -------------------

    double estimated_E = sum_Y / num_samples; // Divide the total sum by the number of samples to get the Monte Carlo estimate.
    return estimated_E; // Return the estimated integral value.
}


#endif //MONTE_CARLO_H
