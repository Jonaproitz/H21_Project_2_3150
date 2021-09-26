
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


void compare_eigenvalues(arma::vec &eigval, arma::vec &lam){
    for (int i = 0; i < lam.size(); i++){
        std::cout << "Eigenvalue " << i + 1 << ":\n    " << eigval(i) << ", " << lam(i)
                  << "\n    " << "Difference: " << fabs(eigval(i) - lam(i)) << "\n\n";
    }
    return;
}


void compare_eigenvectors(arma::mat &eigvec, arma::mat &V){ // Assume normalised eigenvectors
    int N = eigvec.n_cols;
    for (int i = 0; i < N; i++){
        if (arma::approx_equal(V.col(i), -eigvec.col(i), "absdiff", 1e-8)){
            eigvec.col(i) = -eigvec.col(i);
        }
        std::cout << "Eigenvector " << i + 1 << ":\n";
        for (int j = 0; j < N; j++){
            if (eigvec.col(i)(j) < 0){
                std::cout << "    " << eigvec.col(i)(j) << "    ";
            }
            else{
                std::cout << "     " << eigvec.col(i)(j) << "    ";
            }
            if (V.col(i)(j) < 0){
                std::cout << V.col(i)(j);
            }
            else{
                std::cout << " " << V.col(i)(j);
            }
            if (j == 0){
                std::cout << "    Difference: " 
                          << fabs(eigvec.col(i)(j) - V.col(i)(j)) << "\n";
            }
            else{
                std::cout << "                " 
                          << fabs(eigvec.col(i)(j) - V.col(i)(j)) << "\n";
            }
        }
        std::cout << "\n";
    }
    return;
}