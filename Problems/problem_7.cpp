
#include "jacobi_eigensolver.hpp"

int main(){
    int n = 10;
    double L = 1.;

    float dn = L / (n - 1);

    int N = n - 2;

    arma::vec x = arma::linspace(dn, L - dn, N);
}