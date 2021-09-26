

#include "jacobi_eigensolver.hpp"


int main(){
    // Initial constants
    int N = 6;
    double L = 1.;

    double h = L / (N - 1);
    double a = -1 / (h*h);
    double d = 2 / h;

    // Create matrix A of size 6x6
    arma::mat A = create_tridiag_matrix(N, a, d);
    std::cout << "Matrix A =\n";
    A.print();

    // Check analytical solution
    arma::vec lam = arma::vec(N);
    arma::mat V = arma::mat(N, N);
    analytic_solution(lam, V, a, d);


    // Find solution with Jacobi rotation method
    double eps = 1e-8;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1e6, iterations;
    bool converged;
    jacobi_eigensolver(A, eps, eigenvalues, eigenvectors, maxiter, iterations, converged);

    std::cout << "\nEigenvalues calculated by jacobi_eigensolver vs the analytic solution\n";
    compare_eigenvalues(eigenvalues, lam);

    std::cout << "Eigenvectors calculated by jacobi_eigensolver vs the analytic solution\n";
    compare_eigenvectors(eigenvectors, V);
    
    
    std::cout << iterations << "\n";

    std::cout << std::boolalpha;
    std::cout << converged << "\n";

    // End program
    return 0;
}