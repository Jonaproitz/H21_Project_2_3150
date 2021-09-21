#include<iostream>
#include<armadillo>

int main(){
    int N = 6;
    arma::mat A = arma::mat(N, N).fill(0.);
    arma::vec d = arma::vec(6).fill(1.);
    arma::vec a = arma::vec(5). fill(2.);
    A.diag() = d;
    A.diag(1) = a;
    A.diag(-1) = a;
    A.print();
    return 0;
}