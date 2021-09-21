#include<iostream>
#include<armadillo>

int main(){
    int N = 6;
    arma::mat A = arma::mat(N, N).fill(0.);
    A.print();
    return 0;
}