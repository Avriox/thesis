//
// Created by Jakob Goldmann on 05.03.25.
//

#ifndef ONION_MONTE_CARLO_H
#define ONION_MONTE_CARLO_H

#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;
using namespace std;

namespace onion_monte_carlo {
    // Function to calculate binomial coefficient (n choose k)
    double binomial_coefficient(int n, int k) {
        if (k < 0 || k > n) {
            return 0;
        }
        if (k == 0 || k == n) {
            return 1;
        }
        if (k > n / 2) {
            k = n - k;
        }
        double res = 1;
        for (int i = 1; i <= k; ++i) {
            res = res * (n - i + 1) / i;
        }
        return res;
    }

    // Function to calculate double factorial (n!!)
    double double_factorial(int n) {
        if (n < 0) {
            return 1; // Define (-1)!! = 1 and (-3)!! = 1 etc., but for our case we only need for non-negative even numbers and possibly -1!! for j=0 case in moment formula
        }
        if (n == 0 || n == -1) { // Correcting the base case for formula (2j-1)!! when j=0 => -1!! = 1, and for even numbers
            return 1.0;
        }
        double res = 1.0;
        for (int i = n; i > 0; i -= 2) {
            res *= i;
        }
        return res;
    }

    // Function to calculate the (2r)-th moment of a univariate normal distribution N(mu, sigma2)
    double univariate_normal_moment_2r(double mu, double sigma2, int r) {
        double moment = 0.0;
        for (int j = 0; j <= r; ++j) {
            moment += binomial_coefficient(2 * r, 2 * j) * pow(mu, 2 * r - 2 * j) * pow(sigma2, j) * double_factorial(2 * j - 1);
        }
        return moment;
    }


    double analytical_integral_recursive(const vec& mean_vec, const mat& cov_matrix, int r) {
        // TODO dont re-caluclate k. Give k to the function!
        int k = mean_vec.n_elem;

        if (k == 1) {
            return univariate_normal_moment_2r(mean_vec(0), cov_matrix(0, 0), r);
        } else {
            double result_expectation = 0.0;

            // Partition mean vector and covariance matrix
            double mean_beta1 = mean_vec(0);
            vec mean_beta_rest = mean_vec.subvec(1, k - 1);
            double cov_beta1_beta1 = cov_matrix(0, 0);
            rowvec cov_beta1_beta_rest_row = cov_matrix.submat(0, 1, 0, k - 1);
            colvec cov_beta_rest_beta1_col = cov_matrix.submat(1, 0, k - 1, 0);
            mat cov_beta_rest_beta_rest = cov_matrix.submat(1, 1, k - 1, k - 1);

            // Calculate conditional variance and conditional mean function
            mat inv_cov_beta_rest_beta_rest = inv(cov_beta_rest_beta_rest);
            double conditional_variance_beta1 = cov_beta1_beta1 - as_scalar(cov_beta1_beta_rest_row * inv_cov_beta_rest_beta_rest * cov_beta_rest_beta1_col);


            // Define conditional mean as a function of beta_rest:
            auto conditional_mean_beta1_func = [&](const vec& beta_rest) {
                return mean_beta1 + as_scalar(cov_beta1_beta_rest_row * inv_cov_beta_rest_beta_rest * (beta_rest - mean_beta_rest));
            };

            // Calculate E[beta_1^(2r) | beta_{2:k}] which is a function of beta_{2:k}
            auto conditional_expectation_beta1_2r_func = [&](const vec& beta_rest) {
                return univariate_normal_moment_2r(conditional_mean_beta1_func(beta_rest), conditional_variance_beta1, r);
            };

            // Now we need to calculate the expectation of conditional_expectation_beta1_2r_func(beta_{2:k}) * product(beta_{i}^{2r} for i=2 to k)
            // for beta_{2:k} ~ N(mean_beta_rest, cov_beta_rest_beta_rest)
            // Create a new integrand function
            auto integrand_func_recursive = [&](const vec& beta_rest) {
                double product_beta_rest_2r = 1.0;
                for (int i = 0; i < beta_rest.n_elem; ++i) {
                    product_beta_rest_2r *= pow(beta_rest(i), 2.0 * r);
                }
                return conditional_expectation_beta1_2r_func(beta_rest) * product_beta_rest_2r;
            };

            // We cannot analytically integrate this directly in a simple closed form for general r.
            // For comparison, let's use Monte Carlo for the reduced integral over beta_{2:k}
            int num_samples_reduced = 10000000; // Or adjust as needed
            double sum_integrand_values = 0.0;
            mat L_rest = chol(cov_beta_rest_beta_rest, "lower");
            random_device rd{};
            mt19937 gen{rd()};
            normal_distribution<> normal_dist{0, 1};

            for (int i = 0; i < num_samples_reduced; ++i) {
                vec z_rest(k - 1);
                for (int j = 0; j < k - 1; ++j) {
                    z_rest(j) = normal_dist(gen);
                }
                vec beta_rest_sample = mean_beta_rest + L_rest * z_rest;
                sum_integrand_values += integrand_func_recursive(beta_rest_sample);
            }
            result_expectation = sum_integrand_values / num_samples_reduced;

            return result_expectation;
        }
    }
}
#endif //ONION_H
