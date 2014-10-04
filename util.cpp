#include "matrix.h"

double calcResid(Matrix &x, Matrix &A, Matrix &b) 
{
    int N = b.Ni();
//    Matrix resid = b; // this doesn't work!
    Matrix resid = Matrix(N);
    resid = b;
    for(int i=0; i < N; ++i) {
        for(int j=0; j < N; ++j) {
            resid(i) -= A(i,j)*x(j);
        }
    }
    resid *= resid;
    return resid.sum();
}
