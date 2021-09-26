
#include "tridiag_matrix.hpp"

int main(){
    // Setup constants 
    int N = 6;
    double L = 1.;
    double h = L/(N-1);

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

    std::cout << "Eigenvalues calculated by arma::eig_sym vs the analytic solution\n";
    for (int i = 0; i < N; i++){
        std::cout << "Eigenvalue " << i + 1 << ":\n    " << eigval(i) << ", " << lam(i) << "\n";
    }
    
    std::cout << "Eigenvectors calculated by arma::eig_sym vs the analytic solution\n";
    for (int i = 0; i < N; i++){
        std::cout << "Eigenvector " << i + 1 << "\n";
        for (int j = 0; j < N; j++){
            std::cout << "    " << eigvec.col(i)(j) << "    " << V.col(i)(j) << "\n";
        }
    }

    // End program
    return 0;
}