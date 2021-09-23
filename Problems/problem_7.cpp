
#include "jacobi_eigensolver.hpp"

int main(){
    int n = 10;
    double L = 1.;

    float h = L / (n - 1);

    int N = n - 2;

    arma::vec x = arma::linspace(h, L - h, N);


    double a = -1 / (h*h);
    double d = 2 / h;
    arma::mat A = create_tridiag_matrix(N, a, d);
}