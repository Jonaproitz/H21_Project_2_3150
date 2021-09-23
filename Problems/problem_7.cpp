
#include "jacobi_eigensolver.hpp"

int main(){
    // Set constants
    int n = 6;

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

    // Set boundry conditions for v
    arma::vec v_boundry = arma::vec(eigenvalues.size()).fill(0.);

    // Merge solutions with boundries
    x = join_cols(x_0, x, x_n);
    arma::mat V = arma::mat(x.size(), x.size()-1);
    
    arma::vec x_i;
    x_i = x(0);
    V.row(0) = join_rows(x_i, v_boundry.t());
    
    for (int i = 1; i < x.size()-1; i++){
        x_i = {x(i)};
        V.row(i) = join_rows(x_i, eigenvectors.row(i-1));
    }
    x_i = x(x.size()-1);
    V.row(x.size()-1) = join_rows(x_i, v_boundry.t());
    V.print();
    
    // Write solution to binary file
    //arma::mat solution = join_rows(x, V);
    //std::string filename = "solution.bin";
    //solution.save(filename);

    return 0;
}