//
// Created by Jakob Goldmann on 11.03.25.
//

#ifndef ONION_H
#define ONION_H

#include <armadillo>
#include <iostream>
#include <cmath>

using namespace arma;
using namespace std;

namespace onion {


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
        if (r < 0) return 1.0; // for r = 0, moment is 1
        double moment_2r = 0.0;
        for (int j = 0; j <= r; ++j) {
            moment_2r += binomial_coefficient(2 * r, 2 * j) * pow(mu, 2 * r - 2 * j) * double_factorial(2 * j - 1) * pow(sigma2, j);
        }
        return moment_2r;
    }


    double recursive_expected_value(vec mean_vector, mat covariance_matrix, int r) {
        int k_current = mean_vector.n_elem;
        if (k_current == 0) {
            return 1.0; // Product of empty set is 1
        }

        if (k_current == 1) {
            return univariate_normal_moment_2r(mean_vector(0), covariance_matrix(0, 0), r);
        }

        double mu_1 = mean_vector(0);
        vec mu_rest = mean_vector.subvec(1, k_current - 1);
        double Sigma_11 = covariance_matrix(0, 0);
        mat Sigma_1_rest = covariance_matrix.submat(0, 1, 0, k_current - 1);
        mat Sigma_rest_1 = covariance_matrix.submat(1, 0, k_current - 1, 0);
        mat Sigma_rest_rest = covariance_matrix.submat(1, 1, k_current - 1, k_current - 1);

        mat inv_Sigma_rest_rest = inv(Sigma_rest_rest);
        mat cond_cov_mat = Sigma_11 * eye(1,1) - Sigma_1_rest * inv_Sigma_rest_rest * Sigma_rest_1;
        double conditional_variance = cond_cov_mat(0,0);

        // Define a dummy beta_rest (beta_{2:k}) mean to calculate conditional mean as function of beta_{2:k}
        vec beta_rest_val = mu_rest; // For expectation calculation, we use mean of beta_rest

        vec cond_mean_vec = mu_1 * ones<vec>(1) + Sigma_1_rest * inv_Sigma_rest_rest * (beta_rest_val - mu_rest);
        double conditional_mean = cond_mean_vec(0);


        double inner_moment = univariate_normal_moment_2r(conditional_mean, conditional_variance, r);

        if (k_current == 2) {
           vec marginal_mean_beta2 = mu_rest;
           mat marginal_cov_beta2 = Sigma_rest_rest;
           return recursive_expected_value(marginal_mean_beta2, marginal_cov_beta2, r) * inner_moment;
        } else {
            vec marginal_mean_beta_rest = mu_rest;
            mat marginal_cov_beta_rest = Sigma_rest_rest;

            return recursive_expected_value(marginal_mean_beta_rest, marginal_cov_beta_rest, r) * inner_moment;
        }
    }
}
#endif //ONION_H
