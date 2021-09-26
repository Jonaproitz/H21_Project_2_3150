
#include "tridiag_matrix.hpp"

int main(){
    // Setup constants 
    int N = 6;
    double h = 1./(N+1);

    double a = -1/(h*h);
    double d = 2/h;

    // Setup tridiagonal matrix A
    arma::mat A = create_tridiag_matrix(N, a, d);
    std::cout << "Matrix A = \n";
    A.print();

    //Find eigenvalues of A
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A);


    // Find analytical eigenvalues and eigenvectors
    arma::vec lam = arma::vec(N);
    arma::mat V = arma::mat(N, N);
    analytic_solution(lam, V, a, d);

    
    std::cout << "\n\nEigenvalues and eigenvectors calculated by arma::eig_sym vs the analytic solution\n";
    compare_eigen(eigval, eigvec, lam, V);

    // End program
    return 0;
}