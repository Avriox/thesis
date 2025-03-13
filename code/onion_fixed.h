//
// Created by Jakob Goldmann on 12.03.25.
//

#ifndef ONION_FIXED_H
#define ONION_FIXED_H

#include <armadillo>
#include <iostream>
#include <cmath>
#include <map>
#include <vector>
#include <tuple>

using namespace arma;
using namespace std;

namespace onion_fixed {

    // Function to calculate binomial coefficient (n choose k)
    double binomial_coefficient(int n, int k) {
        if (k < 0 || k > n) {
            return 0;
        }
        if (k == 0 || k == n) {
            return 1.0;
        }
        if (k > n / 2) {
            k = n - k;
        }
        double res = 1.0;
        for (int i = 1; i <= k; ++i) {
            res = res * (n - i + 1) / i;
        }
        return res;
    }

    // Function to calculate double factorial
    double double_factorial(int n) {
        if (n < -1) {
            return 0.0; // Or handle as error, double factorial is not defined for negative odd integers
        }
        if (n <= 0) { // Corrected: double_factorial(-1) = 1 and double_factorial(0) = 1 for moment formula consistency
            return 1.0;
        }
        double res = 1.0;
        for (int i = n; i > 0; i -= 2) {
            res *= i;
        }
        return res;
    }


    // Function to calculate univariate normal moment E[X^(2r)] where X ~ N(mu, sigma2)
    double univariate_normal_moment_2r(double mu, double sigma2, int r) {
        if (r <= 0) return 1.0; // for r = 0, moment is 1
        double moment_2r = 0.0;
        for (int j = 0; j <= r; ++j) {
            moment_2r += binomial_coefficient(2 * r, 2 * j) * pow(mu, 2 * r - 2 * j) * double_factorial(2 * j - 1) * pow(sigma2, j);
        }
        return moment_2r;
    }

    // Combined function for conditional mean and variance (for potential minor performance benefit - mostly for single inverse calculation)
    void compute_conditional_mean_var(const arma::mat& mean_vector, const arma::mat& covariance_matrix, double& cond_mean, double& cond_var) {
        // Extract elements from the mean vector and covariance matrix
        double mu_1 = mean_vector(0, 0);  // Mean of the first variable

        // Extract the variance of the first variable
        double sigma_11 = covariance_matrix(0, 0);

        // Extract the covariance between the first variable and the remaining variables
        arma::mat sigma_12 = covariance_matrix.submat(0, 1, 0, covariance_matrix.n_cols - 1);

        // Extract the covariance between the remaining variables and the first variable
        arma::mat sigma_21 = covariance_matrix.submat(1, 0, covariance_matrix.n_rows - 1, 0);

        // Extract the covariance matrix of the remaining variables
        arma::mat sigma_22 = covariance_matrix.submat(1, 1, covariance_matrix.n_rows - 1, covariance_matrix.n_cols - 1);

        // Calculate conditional mean
        // For the case where we're conditioning on the actual values in mean_vector
        cond_mean = mu_1;

        // Calculate conditional variance: Σ₁₁ - Σ₁₂ Σ₂₂⁻¹ Σ₂₁
        arma::mat term = sigma_12 * arma::inv(sigma_22);
        cond_var = sigma_11 - arma::as_scalar(term * sigma_21);
    }

    double get_expected_value(mat mean_vector, mat covariance_matrix, int r) {
        int k = mean_vector.n_elem;

        double cond_mean;
        double cond_var;

        compute_conditional_mean_var(mean_vector, covariance_matrix, cond_mean, cond_var);

        double moment = univariate_normal_moment_2r(cond_mean, cond_var, r);


        return moment;
    }


}

#endif //ONION_FIXED_H
