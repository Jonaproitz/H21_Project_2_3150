
#include "jacobi_eigensolver.hpp"

int main(){
    // Set constants
    int n = 10;

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
    

    // Write solution to binary file
    std::string filename_1 = "solution.bin";
    V.save(filename_1);


    // Call analytic_solution function to find the analytic solution
    arma::vec lam = arma::vec(N);
    arma::mat U = arma::mat(N, N);

    analytic_solution(lam, U, a, d);
    U = join_cols(arma::vec(N).fill(v_boundry(0)).t(), U, arma::vec(N).fill(v_boundry(0)).t());

    // Write analytic solution to binary file
    std::string filename_2 = "analytic_solution.bin";
    U.save(filename_2);


    return 0;
}