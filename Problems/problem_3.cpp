
#include "tridiag_matrix.hpp"

int main(){
    // Setup constants 
    int N = 6;
    double L = 1.;
    double h = L/(N-1);

    double a = -1/h*h;
    double d = 2/h;

    // Setup tridiagonal matrix A
    arma::mat A = create_tridiag_matrix(N, a, d);
    A.print();

    //Find eigenvalues of A
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A);

    eigval.print();
    eigvec.print();


    // Find analytical eigenvalues and eigenvectors
    arma::vec lam = arma::vec(N);
    arma::mat V = arma::mat(N, N);
    analytic_solution(lam, V, a, d);
    
    lam.print();
    V.print();

    // End program
    return 0;
}