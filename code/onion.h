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

// Function to calculate double factorial (for odd double factorial, handle appropriately)
double double_factorial(int n) {
    if (n < 0) {
        return 1.0; // Define (-1)!! = 1 for the formula to work when j=0
    }
    if (n == 0) {
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
    double moment_2r = 0.0;
    for (int j = 0; j <= r; ++j) {
        moment_2r += binomial_coefficient(2 * r, 2 * j) * pow(mu, 2 * r - 2 * j) * double_factorial(2 * j - 1) * pow(sigma2, j);
    }
    return moment_2r;
}

// Recursive function to calculate E[prod_{i=1}^k beta_i^(2r)]
double recursive_expected_value(vec mean_vector, mat covariance_matrix, int r) {
    // Base case: 1-dimensional distribution
    if (mean_vector.n_elem == 1) {
        return univariate_normal_moment_2r(mean_vector(0), covariance_matrix(0, 0), r);
    }

    int k_current = mean_vector.n_elem; // Current dimension

    // 1. Partition parameters
    double mu_1 = mean_vector(0);
    vec mu_rest = mean_vector.subvec(1, k_current - 1);
    double Sigma_11 = covariance_matrix(0, 0);
    mat Sigma_1_rest, Sigma_rest_1, Sigma_rest_rest;

    if (k_current > 1) {
        Sigma_1_rest = covariance_matrix.submat(0, 1, 0, k_current - 1);
        Sigma_rest_1 = covariance_matrix.submat(1, 0, k_current - 1, 0);
        Sigma_rest_rest = covariance_matrix.submat(1, 1, k_current - 1, k_current - 1);
    }


    // 2. Calculate conditional variance
    double conditional_variance = Sigma_11; // Initialize in case k_current = 1
    if (k_current > 1) {
        mat inv_Sigma_rest_rest = inv(Sigma_rest_rest);
        mat conditional_variance_mat = Sigma_11 * eye(1,1) - Sigma_1_rest * inv_Sigma_rest_rest * Sigma_rest_1;
        conditional_variance = conditional_variance_mat(0,0);
    }


    double expected_value = 0.0;

    // 3. Sum over j from 0 to r
    for (int j = 0; j <= r; ++j) {
        double v_j = pow(conditional_variance, j) * double_factorial(2 * j - 1);
        int power_of_conditional_mean = 2 * r - 2 * j;
        double binomial_coeff = binomial_coefficient(2 * r, 2 * j);

        double mean_conditional_mean = mu_1; // Approximation
        double mean_conditional_mean_term = pow(mean_conditional_mean, power_of_conditional_mean);

        double remaining_expectation = 1.0;
        if (k_current > 1) {
            remaining_expectation = recursive_expected_value(mu_rest, Sigma_rest_rest, r); // Recursive call
        }

        expected_value += binomial_coeff * mean_conditional_mean_term * v_j * remaining_expectation;
    }
    return expected_value;
}

#endif //ONION_H
