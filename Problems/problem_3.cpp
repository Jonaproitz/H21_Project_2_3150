#include<iostream>
#include<armadillo>

int main(){
    int N = 6;
    arma::mat A = arma::mat(N, N).fill(0.);
    arma::vec d = arma::vec(6).fill(1.);
    A.diag() = d;
    A.print();
    return 0;
}