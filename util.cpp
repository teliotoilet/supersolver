#include "matrix.h"

double calcResid(Matrix &x, Matrix &A, Matrix &b) 
{
    int N = b.Ni();
    for(int i=0; i < N; ++i) {
        for(int j=0; j < N; ++j) {
            b(i) -= A(i,j)*x(j);
        }
    }
    b *= b;
    return b.sum();
}
