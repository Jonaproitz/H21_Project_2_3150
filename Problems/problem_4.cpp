
#include "max_offdiag_symmetric.hpp"

int main(){
    // Call max_offdiagonal_symmetric function
    arma::mat A = {
        {1., 0., 0., 0.5},
        {0., 1., -0.7, 0.},
        {0., -0.7, 1., 0.},
        {0.5, 0., 0., 1.}
        };
    int k;
    int l;
    double a = max_offdiag_symmetric(A, k ,l);
    std::cout << a << "\n";
    
    // End program
    return 0;
}