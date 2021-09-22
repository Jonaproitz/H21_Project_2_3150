

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

    
    // Find original indices of largest off-diagonal matrix element and abs(a_max)
    

    // Loop Jacobi_rotate until abs(a_max)<eps


        // Stops if the number of iterations reaches "maxiter"


    // Store number of iterations in integer "iterations"


    // Set bool reference "converged" to true


    // Store eigenvalues and eigenvectors


    // End function with no return
    return;
}


int main(){

    // End program
    return 0;
}