#include<iostream>
#include<armadillo>

int main(){
    // Setup constants 
    int N = 6;
    double L = 1.;
    double h = L/(N-1);

    double a = -1/h*h;
    double d = 2/h;

    // Setup tridiagonal matrix A
    arma::mat A = arma::mat(N, N).fill(0.);
    arma::vec b = arma::vec(N).fill(d);
    arma::vec ac = arma::vec(N-1). fill(a);
    A.diag() = b;
    A.diag(1) = ac;
    A.diag(-1) = ac;
    A.print();

    // End program
    return 0;
}