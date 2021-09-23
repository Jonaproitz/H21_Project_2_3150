
#include "jacobi_eigensolver.hpp"

int main(){
    // Set constants
    int n = 6;

    float h = 1. / n;

    int N = n - 1;

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


    // Set boundry conditions for v
    arma::vec v_boundry = {0.};//arma::vec(eigenvalues.size()).fill(0.);

    // Merge solutions with boundries
    x = join_cols(x_0, x, x_n);
    arma::mat V = arma::mat(x.size(), x.size()-1);

    V.col(0) = x;
    for (int i = 1; i < x.size()-1; i++){
        V.col(i) = join_cols(v_boundry, eigenvectors.col(i-1), v_boundry);
    }
    V.print();
    

    // Write solution to binary file
    std::string filename = "solution.bin";
    V.save(filename);

    return 0;
}