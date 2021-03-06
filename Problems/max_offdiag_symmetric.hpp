
#include "tridiag_matrix.hpp"

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l){
    int n = A.n_rows;
    double a_max;

    for (int i = 0; i < n; i++){
        for (int j = i+1; j < n; j++){
            double val = fabs(A(i, j));
            if (val > a_max){
                a_max = val; k = i; l = j;
            }
        }
    }
    
    // Return matrix element
    return a_max;
}