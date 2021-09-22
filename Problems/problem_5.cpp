

#include "max_offdiag_symmetric.hpp"


void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    // Necessary computations before rotation
    double t, c, s;
    if (A(k, l) != 0){
        double tau = (A(l ,l) - A(k, k)) / (2*A(k, l));
        if (tau >= 0){
            t = 1 / (tau + sqrt(1 + tau*tau));
        }
        else{
            t = -1 / (-tau + sqrt(1 + tau*tau));
        }

        double c = 1 / sqrt(1 + t*t);
        double s = c*t;
    }
    else{
        t = 0;
        c = 1;
        s = 0;
    }


    // Perform a single jacobi rotation
    double a_kk = A(k ,k), a_ll = A(l, l), a_kl = A(k, l);
    A(k, k) = a_kk*c*c - 2*a_kl*c*s + a_ll*s*s;
    A(l ,l) = a_ll*c*c + 2*a_kl*c*s + a_kk*s*s;
    A(k, l) = 0;
    A(l, k) = 0;

    double a_ik, a_il;
    for (int i = 0; i < A.n_rows; i++){
        if (i != k && i != l){
            a_ik = A(i, k);
            a_il = A(i, l);

            A(i, k) = a_ik * c - a_il*s;
            A(k, i) = A(i, k);

            A(i, l) = a_il*c - a_ik*s;
            A(l, i) = A(i, l);
        }
    }
    double r_ik, r_il;
    for (int i = 0; i < A.n_rows; i++){
        r_ik = R(i, k);
        r_il = R(i, l);
        R(i, k) = r_ik*c - r_il*s;
        R(i, l) = r_il*c - r_ik*s;
    }

    // End funtion with no return
    return;
}


void jacobi_eigensolver(const arma::mat& A, double eps, arma::vec& eigenvalues, 
                        arma::mat& eigenvectors, const int maxiter, int& iterations, 
                        bool& converged){
    // Create necessary inputs for jacobi_rotate
    arma::mat A_m = A;
    arma::mat R;
    R.eye(size(A));

    
    // Find original indices of largest off-diagonal matrix element and abs(a_max)
    int k, l;
    double a_max;
    a_max = max_offdiag_symmetric(A, k, l);


    // Loop Jacobi_rotate until abs(a_max)<eps
    int i = 0;
    while (a_max > eps && i < maxiter){
        jacobi_rotate(A_m, R, k, l);
        a_max = max_offdiag_symmetric(A_m, k, l);

        // Stops if the number of iterations reaches "maxiter"
        i++;
    }


    // Store number of iterations in integer "iterations"
    iterations = i + 1;


    // Set bool reference "converged" to true
    converged = true;


    // Store eigenvalues and eigenvectors
    eigenvalues = A_m.diag(0);
    eigenvectors = R;

    // End function with no return
    return;
}


int main(){
    // Initial constants
    int N = 6;
    double L = 1.;

    double h = L / (N - 1);
    double a = -1 / (h*h);
    double d = 2 / h;

    // Create matrix A of size 6x6
    arma::mat A = create_tridiag_matrix(N, a, d);
    
    // Check analytical solution
    arma::vec lam = arma::vec(N);
    arma::mat V = arma::mat(N, N);
    analytic_solution(lam, V, a, d);

    lam.print();
    V.print();


    // Find solution with Jacobi rotation method
    double eps = 1e-8;
    arma::vec eigenvalues;
    arma::mat eigenvectors;
    int maxiter = 1e6;
    int iterations;
    bool converged;


    // End program
    return 0;
}