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
    arma::vec ac = arma::vec(N-1).fill(a);
    A.diag() = b;
    A.diag(1) = ac;
    A.diag(-1) = ac;
    A.print();

    //Find eigenvalues of A
    arma::vec eigval;
    arma::mat eigvec;

    arma::eig_sym(eigval, eigvec, A);

    eigval.print();
    eigvec.print();


    // Find analytical eigenvalues and eigenvectors
    double pi = arma::datum::pi;
    arma::vec lam = arma::vec(N);
    arma::mat V = arma::mat(N, N);
    for (int i = 1; i <= N; i++){
        lam(i-1) = d + 2*a * cos(i * pi / (N + 1));
        for (int j = 1; j <= N; j++){
            V(i-1, j-1) = sin(j * i * pi / (N + 1));
        }
    }
    arma::mat V_norm = arma::normalise(-V);
    lam.print();
    V_norm.print();

    // End program
    return 0;
}