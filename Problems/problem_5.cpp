

#include "max_offdiag_symmetric.hpp"


void jacobi_rotate(arma::mat& A, arma::mat& R, int k, int l){
    // Necessary computations before rotation
    double t, c, s;

    if (A(k, l) != 0){
        double tau = (A(l ,l) - A(k, k)) / (2*A(k, l));
        if (tau > 0){
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