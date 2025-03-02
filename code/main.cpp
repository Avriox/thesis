#include <armadillo>
#include <iostream>

int main() {

    // -------------------------
    // *** Initial Variables ***
    // -------------------------
    
    // Set matrix dimensions (n rows × k columns) for X_k
    const int n = 5;  
    const int k = 3; 
    
    // Verify n > k
    if (n <= k) {
        std::cerr << "Error: n must be greater than k" << std::endl;
        return 1;
    }

    // Create random matrix with values between 0 and 1
    arma::arma_rng::set_seed(123);
    arma::mat X_k = arma::randu<arma::mat>(n, k);

    // For now tau is some positive real number
    double tau = 0.2;

    // For now, A_k is the identity matrix of size k
    arma::mat A_k = arma::eye<arma::mat>(k, k);


    // -------------------
    // *** Computation ***
    // -------------------
    
    // Ck = Transpose[Xk].Xk + (1 / τ) * Ak

    arma::mat test = X_k * X_k.t();
    test.print();
    
    // arma::mat C_k = X_k.t() * X_k + (1.0/tau) * A_k;

    // C_k.print();
    
    return 0;
}