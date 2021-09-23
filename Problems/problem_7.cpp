
#include "jacobi_eigensolver.hpp"

int main(){
    int n = 10;

    float h = 1. / (n - 1);

    int N = n - 2;

    arma::vec x = arma::linspace(h, 1. - h, N);


    double a = -1 / (h*h);
    double d = 2 / h;
    arma::mat A = create_tridiag_matrix(N, a, d);


    double eps = 1e-8;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1e6, iterations;
    bool converged;
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    return 0;
}