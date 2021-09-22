

#include "jacobi_eigensolver.hpp"


int main(int argc, char* argv[]){
    // Add command line argument N
    if (argc != 1){
        std::cout << "Please input 1 value";
        return 1;
    }
    int N = atoi(argv[1]);

    // Calculate variables for diagonals
    double L = 1.;
    double h = L / (N - 1);
    
    double a = -1 / (h*h);
    double d = 2 / h;

    // Decleare necessary vairables
    double eps = 1e-8;
    arma::vec eigenvalues;
    arma::mat A = create_tridiag_matrix(N, a, d), eigenvectors;
    int maxiter = 1e6, iterations;
    bool converged;

    // Call eigensolver
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);


    return 0;
}