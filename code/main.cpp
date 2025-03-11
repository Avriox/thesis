#include <armadillo>
#include <iostream>
#include "monte_carlo.h"
#include "onion.h"
#include "onion_monte_carlo.h"
#include "timer.h"

int main() {

    // -------------------------
    // *** Initial Variables ***
    // -------------------------
    
    // Set matrix dimensions (n rows × k columns) for X_k
    const int n = 15;
    const int k = 10;
    
    // Verify n > k
    if (n <= k) {
        std::cerr << "Error: n must be greater than k" << std::endl;
        return 1;
    }

    // Create random matrix with values between 0 and 1
    arma::arma_rng::set_seed(123);
    arma::mat X_k = arma::randu(n, k);

    // For now tau is some positive real number
    double tau = 0.2;

    // For now, A_k is the identity matrix of size k
    arma::mat A_k = arma::eye<arma::mat>(k, k);

    // y is a "data vector of size n". Random for now
    arma::vec y = arma::randu(n);

    double sigma_squared = 0.56;
    double r = 6;

    // -------------------
    // *** Computation ***
    // -------------------
    
    // Ck = Transpose[Xk].Xk + (1 / τ) * Ak
    arma::mat C_k = X_k.t() * X_k + (1.0/tau) * A_k;

    // β₍ₖ₎ = C₍ₖ₎⁻¹ X₍ₖ₎ᵀ y₍ₙ₎
    // TODO use a more efficient decomposition to get inv and det of C_k.
    arma::mat beta_tk =  C_k.i() * X_k.t() * y;

    // R₍ₖ₎=y₍ₙ₎ᵀ(I₍ₙ₎ - X₍ₖ₎ C₍ₖ₎⁻¹ X₍ₖ₎ᵀ) y₍ₙ₎
    // TODO re-use the inverse of C_k! Reuse X_k transpose!
    arma::mat R_k = y.t() * (arma::eye<arma::mat>(n, n) - X_k * C_k.i() * X_k.t()) * y;

    // // --- Monte Carlo Method ---
    // double estimated_integral_mc = TIME(
    //     "Monte Carlo Integration",
    //     solve_integral_monte_carlo_cpp(C_k, X_k, y, sigma_squared, r, 10000000)
    //     );
    // std::cout << "Estimated integral value using Monte Carlo: " << estimated_integral_mc << std::endl;
    //
    // mat Sigma = sigma_squared * C_k.i();
    //
    // // --- Recursive Method ---
    // double estimated_integral_analytical_hybrid = TIME(
    //     "Analytical / Montecarlo Hybrid",
    //     onion_monte_carlo::analytical_integral_recursive(beta_tk, Sigma, r)
    // );
    // cout << "Estimated integral value using Recursive (Hybrid MC) approach: " << estimated_integral_analytical_hybrid << endl;


    mat Sigma = sigma_squared * C_k.i();

    // --- Onion Method ---
    double estimated_integral_analytical = TIME(
        "Analytical ONION",
        recursive_expected_value(beta_tk, Sigma, r)
    );
    cout << "Estimated integral value using Recursive onion approach: " << estimated_integral_analytical << endl;


    timer_report();

    return 0;
}
