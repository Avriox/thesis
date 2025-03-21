#include <armadillo>
#include <iostream>
#include "onion_fixed.h"
#include "timer.h"
#include "isserlis_theorem.h"

int main() {

    // -------------------------
    // *** Initial Variables ***
    // -------------------------
    
    // Set matrix dimensions (n rows × k columns) for X_k
    const int n = 1000;   //zehn tausende aber mal 1000 bis 10000
    const int k = 10;    // 50 - 500 
    
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
    double r = 2;

    // -------------------
    // *** Computation ***
    // -------------------
    
    // Ck = Transpose[Xk].Xk + (1 / τ) * Ak
    arma::mat C_k = X_k.t() * X_k + (1.0/tau) * A_k;

    // β₍ₖ₎ = C₍ₖ₎⁻¹ X₍ₖ₎ᵀ y₍ₙ₎     /* Renamed to mu for mean
    // TODO use a more efficient decomposition to get inv and det of C_k.
    arma::mat mu =  C_k.i() * X_k.t() * y;

    // R₍ₖ₎=y₍ₙ₎ᵀ(I₍ₙ₎ - X₍ₖ₎ C₍ₖ₎⁻¹ X₍ₖ₎ᵀ) y₍ₙ₎
    // TODO re-use the inverse of C_k! Reuse X_k transpose!
    arma::mat R_k = y.t() * (arma::eye<arma::mat>(n, n) - X_k * C_k.i() * X_k.t()) * y;

    mat Sigma = sigma_squared * C_k.i();

    double expected_value_mombf = TIME(
        "mombf",
        mombf::get_expected_value(mu, Sigma, r)
    );
    cout << "Expected value using mombf approach: " << expected_value_mombf << endl;
    
    // double expected_value_mombfarma = TIME(
    //     "mombf arma",
    //     mombf_arma::get_expected_value(mu, Sigma, r)
    // );
    // cout << "Expected value using mombf approach: " << expected_value_mombfarma << endl;
    
    
    // --- Onion Method  ---
    // double expected_value_isserils = TIME(
    //     "Onion",
    //     isserlis::get_expected_value(mu, Sigma, r)
    // );
    // cout << "Expected value using onion approach: " << expected_value_isserils << endl;
    
    
    // --- Onion Method  ---
    // double expected_value = TIME(
    //     "Onion",
    //     onion_fixed::get_expected_value(mu, Sigma, r)
    // );
    // cout << "Expected value using onion approach: " << expected_value << endl;





    timer_report();

    return 0;
}
