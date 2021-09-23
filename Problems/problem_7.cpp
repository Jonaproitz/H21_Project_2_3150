
#include "jacobi_eigensolver.hpp"

int main(){
    // Set constants
    int n = 10;

    float h = 1. / (n - 1);

    int N = n - 2;

    // Create x^ vector and boundry conditions
    arma::vec x = arma::linspace(h, 1. - h, N);
    arma::vec x_0 = {0.};
    arma::vec x_n = {1.};

    // Calculate more constants
    double a = -1 / (h*h);
    double d = 2 / h;
    arma::mat A = create_tridiag_matrix(N, a, d);

    // Solve DE with Jacobi rotation algorithm
    double eps = 1e-8;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1e6, iterations;
    bool converged;
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    // Set boundryconditions for v
    arma::vec v_boundry = {0.};

    


    return 0;
}