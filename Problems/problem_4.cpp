
#include <armadillo>
#include <iostream>

double max_offdiag_symmetric(const arma::mat& A, int& k, int &l){
    int n = A.col(0).n_elem - 1;
    double a_max;

    for (int i = 1; i < n; i++){
        for (int j = i+1; j < n; j++){
            double val = fabs(A(i, j));
            if (val > a_max){
                a_max = val; k = i; l = j;
            }
        }
    }
    
    // Return matrix element
    return A(k, l);
}

int main(){
    // Call max_offdiagonal_symmetric function
    arma::mat A = {
        {1., 0., 0., 0.5},
        {0., 1., -0.7, 0.},
        {0., -0.7, 1., 0.},
        {0.5, 0., 0., 1.}
        };
    int k = 1;
    int l = 1;
    double a = max_offdiag_symmetric(A, k ,l);
    std::cout << a << "\n";
    
    // End program
    return 0;
}