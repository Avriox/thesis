//
// Created by jakob on 3/14/25.
//

#ifndef ISSERLIS_THEOREM_H
#define ISSERLIS_THEOREM_H
#include <armadillo_bits/typedef_mat.hpp>
#include <armadillo_bits/typedef_mat.hpp>

namespace isserlis {

    double expected_product(arma::ivec indices, const arma::vec& mean, const arma::mat& cov) {
        if (indices.n_elem == 0) {
            return 1.0;
        }

        if (indices.n_elem == 2) {
            int index1 = indices(0);
            int index2 = indices(1);
            return cov(index1, index2) + mean(index1) * mean(index2);
        }

        double result = 0.0;
        int first_index = indices(0);
        arma::ivec remaining_indices = indices.subvec(1, indices.n_elem - 1);

        for (int i = 0; i < remaining_indices.n_elem; ++i) {
            int second_index = remaining_indices(i);
            double expected_pair = cov(first_index, second_index) + mean(first_index) * mean(second_index);

            arma::ivec next_indices;
            if (i == 0) {
                next_indices = remaining_indices.subvec(1, remaining_indices.n_elem - 1);
            } else if (i == remaining_indices.n_elem - 1) {
                next_indices = remaining_indices.subvec(0, remaining_indices.n_elem - 2);
            } else {
                arma::ivec part1 = remaining_indices.subvec(0, i - 1);
                arma::ivec part2 = remaining_indices.subvec(i + 1, remaining_indices.n_elem - 1);
                next_indices = arma::join_cols(part1, part2);
            }

            result += expected_pair * expected_product(next_indices, mean, cov);
        }

        return result;
    }

    double get_expected_value(arma::vec mean, arma::mat cov, int r) {
        int k = mean.n_rows;

        arma::ivec indices = arma::ivec(k * 2 * r);
        int current_index = 0;
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < 2 * r; ++j) {
                indices(current_index++) = i;
            }
        }


        return expected_product(indices, mean, cov);
    }

    
}

namespace mombf_arma {
    // Helper function to get the next n-tuple in lexicographical order
    int GetNextTuple(arma::ivec& tuple, int n, int base) {
        int j = 0;
        while (j < n && tuple(j) == base - 1) {
            tuple(j) = 0;
            j++;
        }
        if (j < n)
            tuple(j)++;
        return (j < n);
    }

    // Helper function to calculate the binomial coefficient (limited to power 2 or 4 as in the original code)
    int BinomialCoefficient(int power, int nu) {
        int ans;
        if (power == 2) {
            ans = 1 + nu % 2;
        } else if (power == 4) {
            if (nu == 0 || nu == 4)
                ans = 1;
            if (nu == 1 || nu == 3)
                ans = 4;
            if (nu == 2)
                ans = 6;
        } else
            ans = 0; /* should never happen based on original code's usage */
        return ans;
    }

    // Helper function to calculate one_plus_kappa (specific to t-distribution)
    double one_plus_kappa(double dof, int r) {
        double product = 1.0;
        int i;

        if (r == 0)
            return 1.0;

        for (i = 1; i <= r; i++)
            product *= 0.5 * dof - i;
        return pow(0.5 * dof - 1.0, r) / product;
    }

    // Helper function to calculate the logarithm of the factorial
    double lfact(int n) {
        if (n < 0)
            return 0.0; // Or handle error appropriately
        if (n <= 1)
            return 0.0;
        double sum = 0.0;
        for (int i = 2; i <= n; i++)
            sum += log(static_cast<double>(i));
        return sum;
    }

    // Function adapted from the provided C code to use Armadillo vectors and matrices
    double mvtexpect_arma(const arma::vec& mu_arma, const arma::mat& sigma_arma, int n, int power, double dof) {
        double product = 1.0;
        double sum = 0.0;
        double covariance_term, mean_term, temp;
        int index_sum;

        int s = power * n;
        int half_s = s / 2;
        int j, k, r, half_power;

        arma::ivec nu_index = arma::zeros<arma::ivec>(n); /* indices nu_0 through nu_{n-1} */

        for (r = 0; r <= half_s; r++) {
            for (j = 0; j < n; j++)
                nu_index(j) = 0;

            do {
                product = 1.0;

                index_sum = 0;
                for (j = 0; j < n; j++)
                    index_sum += nu_index(j);
                if (index_sum % 2)
                    product *= -1;

                for (j = 0; j < n; j++)
                    product *= BinomialCoefficient(power, nu_index(j));

                if (dof > 0.0)
                    product *= one_plus_kappa(dof, r);
                /* else normal case and multiplicative term is 1 */

                covariance_term = 0.0;
                half_power = power / 2;
                for (j = 0; j < n; j++) {
                    temp = 0.0; /* double_sum accumulates the product sigma*h */
                    for (k = 0; k < n; k++) {
                        /* sigma is stored by columns in the original code, so sigma_{ij} = sigma[i + j*n] */
                        /* Assuming sigma_arma is row-major, we adjust the access accordingly */
                        temp += sigma_arma(j, k) * (half_power - nu_index(k));
                    }
                    covariance_term += (half_power - nu_index(j)) * temp;
                }
                product *= pow(0.5 * covariance_term, r);

                mean_term = 0.0;
                for (j = 0; j < n; j++)
                    mean_term += (half_power - nu_index(j)) * mu_arma(j); // Adjusted index

                product *= pow(mean_term, s - 2 * r);
                product /= exp(lfact(r) + lfact(s - 2 * r));

                sum += product;

            } while (GetNextTuple(nu_index, n, power + 1));
        }
        return sum;
    }

    double get_expected_value(const arma::vec& mean, const arma::mat& cov, int r) {
        int k = mean.n_rows;
        double dof = -1.0; // For normal distribution
        return mvtexpect_arma(mean, cov, k, r, dof);
    }
}

namespace mombf {
    // Helper function to allocate an integer vector (mimicking ivector from the C code)
    int* ivector(int nl, int nh) {
        int* v = new int[nh - nl + 2]; // Adjusted size to match C indexing
        return v - nl;
    }

    // Helper function to free an integer vector (mimicking free_ivector from the C code)
    void free_ivector(int* v, int nl, int nh) {
        delete(v + nl);
    }

    // Helper function to get the next n-tuple in lexicographical order
    int GetNextTuple(int* tuple, int n, int base) {
        int j = 0;
        while (j < n && tuple[j] == base - 1) {
            tuple[j] = 0;
            j++;
        }
        if (j < n)
            tuple[j]++;
        return (j < n);
    }

    // Helper function to calculate the binomial coefficient (limited to power 2 or 4 as in the original code)
    int BinomialCoefficient(int power, int nu) {
        int ans;
        if (power == 2) {
            ans = 1 + nu % 2;
        } else if (power == 4) {
            if (nu == 0 || nu == 4)
                ans = 1;
            if (nu == 1 || nu == 3)
                ans = 4;
            if (nu == 2)
                ans = 6;
        } else
            ans = 0; /* should never happen based on original code's usage */
        return ans;
    }

    // Helper function to calculate one_plus_kappa (specific to t-distribution)
    double one_plus_kappa(double dof, int r) {
        double product = 1.0;
        int i;

        if (r == 0)
            return 1.0;

        for (i = 1; i <= r; i++)
            product *= 0.5 * dof - i;
        return pow(0.5 * dof - 1.0, r) / product;
    }

    // Helper function to calculate the logarithm of the factorial
    double lfact(int n) {
        if (n < 0)
            return 0.0; // Or handle error appropriately
        if (n <= 1)
            return 0.0;
        double sum = 0.0;
        for (int i = 2; i <= n; i++)
            sum += log(static_cast<double>(i));
        return sum;
    }

    // Function adapted from the provided C code
    double mvtexpect_adapted(double* mu, double** sigma, int n, int power, double dof) {
        double product = 1.0;
        double sum = 0.0;
        double covariance_term, mean_term, temp;
        int index_sum;

        int s = power * n;
        int half_s = s / 2;
        int j, k, r, half_power;

        int* nu_index; /* indices nu_0 through nu_{n-1} */
        nu_index = ivector(0, n - 1); // Adjusted index range

        for (r = 0; r <= half_s; r++) {
            for (j = 0; j < n; j++)
                nu_index[j] = 0;

            do {
                product = 1.0;

                index_sum = 0;
                for (j = 0; j < n; j++)
                    index_sum += nu_index[j];
                if (index_sum % 2)
                    product *= -1;

                for (j = 0; j < n; j++)
                    product *= BinomialCoefficient(power, nu_index[j]);

                if (dof > 0.0)
                    product *= one_plus_kappa(dof, r);
                /* else normal case and multiplicative term is 1 */

                covariance_term = 0.0;
                half_power = power / 2;
                for (j = 0; j < n; j++) {
                    temp = 0.0; /* double_sum accumulates the product sigma*h */
                    for (k = 0; k < n; k++) {
                        /* sigma is stored by columns in the original code, assuming row-major here as per typical C++ array */
                        temp += sigma[j][k] * (half_power - nu_index[k]);
                    }
                    covariance_term += (half_power - nu_index[j]) * temp;
                }
                product *= pow(0.5 * covariance_term, r);

                mean_term = 0.0;
                for (j = 0; j < n; j++)
                    mean_term += (half_power - nu_index[j]) * mu[j]; // Adjusted index

                product *= pow(mean_term, s - 2 * r);
                product /= exp(lfact(r) + lfact(s - 2 * r));

                sum += product;

            } while (GetNextTuple(nu_index, n, power + 1));
        }
        free_ivector(nu_index, 0, n - 1); // Adjusted index range
        return sum;
    }

    double get_expected_value(const arma::vec& mean_arma, arma::mat cov_arma, int r) {
        int k = mean_arma.n_rows;
       

        // Convert arma::vec to double*
        double* mu = new double[k];
        for (int i = 0; i < k; ++i) {
            mu[i] = mean_arma(i);
        }

        // Convert arma::mat to double**
        double** sigma = new double*[k];
        for (int i = 0; i < k; ++i) {
            sigma[i] = new double[k];
            for (int j = 0; j < k; ++j) {
                sigma[i][j] = cov_arma(i, j);
            }
        }

        double dof = -1.0; // For normal distribution
        double result = mvtexpect_adapted(mu, sigma, k, r, dof);

        // Clean up allocated memory
        delete mu;
        for (int i = 0; i < k; ++i) {
            delete sigma[i];
        }
        delete sigma;

        return result;
    }
}
#endif //ISSERLIS_THEOREM_H
