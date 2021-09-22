
#include <armadillo>
#include <iostream>


arma::mat create_tridiag_matrix(int N, double a, double d){
    // Create tridiagonal matrix A
    arma::mat A = arma::mat(N, N).fill(0.);
    arma::vec b = arma::vec(N).fill(d);
    arma::vec ac = arma::vec(N-1).fill(a);
    A.diag(0) = b;
    A.diag(1) = ac;
    A.diag(-1) = ac;

    // Return tridiagonal matrix A
    return A;
}


void analytic_solution(arma::vec &lam, arma::mat &V, double a, double d){
    double pi = arma::datum::pi;
    int N = lam.n_elem;
    for (int i = 1; i <= N; i++){
        lam(i-1) = d + 2*a * cos(i * pi / (N + 1));
        for (int j = 1; j <= N; j++){
            V(i-1, j-1) = sin(j * i * pi / (N + 1));
        }
    }
    V = arma::normalise(-V);

    // End function
    return;
}