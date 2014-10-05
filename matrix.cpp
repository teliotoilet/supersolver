#include "main.h"
#include "matrix.h"

#ifdef DEBUG
#include <cassert>
#endif

Matrix::Matrix(int ni, int nj, int nk)
    : _ni(ni), _nj(nj), _nk(nk)
{
    _n = ni*nj*nk;
    _A.resize(_n);
    if(verbose) std::cout << "Constructed matrix with dimensions (" 
        << ni << "," << nj << "," << nk << ")\n";
}
Matrix::~Matrix()
{
    if(verbose) std::cout << "destroyed matrix with dimensions ("
        << _ni << "," << _nj << "," << _nk << ")\n";
}

double& Matrix::operator() (const int i)
{
#ifdef DEBUG
    assert(i >= 0 && i < _ni);
#endif
    return _A[i];
}
double& Matrix::operator() (const int i, const int j)
{
#ifdef DEBUG
    assert(i >= 0 && i < _ni);
    assert(j >= 0 && j < _nj);
#endif
    return _A[i+_ni*j];
}
double& Matrix::operator() (const int i, const int j, const int k)
{
#ifdef DEBUG
    assert(i >= 0 && i < _ni);
    assert(j >= 0 && j < _nj);
    assert(k >= 0 && k < _nk);
#endif
    return _A[i + _ni*j + _ni*_nj*k];
}

std::ostream& operator<< (std::ostream &out, Matrix &M)
{
    out << "[";
    for(int i=0; i < M._ni; ++i) {
        if(i==0) out << "[";
        else out << ",";

        for(int j=0; j < M._nj; ++j) {
            if(j==0) out << "[";
            else out << ",";

            for(int k=0; k < M._nk; ++k) {
                if(k==0) out << "[";
                else out << ",";

                out << M(i,j,k);

            }
            out << "]\n";
        }
        out << "]";
    }
    out << "]";
}

double Matrix::sum() {
    double d=0.0;
    for(int idx=0; idx < _n; ++idx) {
        d += _A[idx];
    }
    return d;
}

Matrix& Matrix::operator= (const Matrix &m)
{
#ifdef DEBUG
    assert(_ni == m._ni);
    assert(_nj == m._nj);
    assert(_nk == m._nk);
#endif
    for(int idx=0; idx < _n; ++idx) {
        _A[idx] = m._A[idx];
    }
    //if(verbose) std::cout << "Overloaded assignment (const) with dimensions (" 
    //    << _ni << "," << _nj << "," << _nk << ")\n";
    return *this;
}
Matrix& Matrix::operator= (const double d)
{
    for(int idx=0; idx < _n; ++idx) {
        _A[idx] = d;
    }
    return *this;
}
Matrix& Matrix::operator= (const double d[])
{
    for(int idx=0; idx < _n; ++idx) {
        _A[idx] = d[idx];
    }
    if(verbose) std::cout << "Overloaded assignment (array) with dimensions (" 
        << _ni << "," << _nj << "," << _nk << ")\n";
    return *this;
}

//Matrix operator+(const Matrix &a1, const Matrix &a2)
//{
//#ifdef DEBUG
//    assert(a1._ni==a2._ni);
//    assert(a1._nj==a2._nj);
//    assert(a1._nk==a2._nk);
//#endif
//    int ni = a1._ni;
//    int nj = a1._nj;
//    int nk = a1._nk;
//    Matrix A = Matrix(ni,nj,nk);
//    for(int i=0; i < ni; ++i) {
//        for(int j=0; j < nj; ++j) {
//            for(int k=0; k < nk; ++k) {
//                A(i,j,k) = a1._A[i][j][k] + a2._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
//Matrix operator+(const Matrix &m, const double d)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = m._A[i][j][k] + d;
//            }
//        }
//    }
//    return A;
//}
//Matrix operator+(const double d, const Matrix &m)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = d + m._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
//Matrix& Matrix::operator+=(const Matrix &m)
//{
//    for(int i=0; i < _ni; ++i) {
//        for(int j=0; j < _nj; ++j) {
//            for(int k=0; k < _nk; ++k) {
//                _A[i][j][k] += m._A[i][j][k];
//            }
//        }
//    }
//    return *this;
//}
//Matrix& Matrix::operator+=(const double d)
//{
//    for(int i=0; i < _ni; ++i) {
//        for(int j=0; j < _nj; ++j) {
//            for(int k=0; k < _nk; ++k) {
//                _A[i][j][k] += d;
//            }
//        }
//    }
//    return *this;
//}
//Matrix& Matrix::operator-=(const Matrix &m)
//{
//    for(int i=0; i < _ni; ++i) {
//        for(int j=0; j < _nj; ++j) {
//            for(int k=0; k < _nk; ++k) {
//                _A[i][j][k] -= m._A[i][j][k];
//            }
//        }
//    }
//    return *this;
//}
//Matrix& Matrix::operator-=(const double d)
//{
//    for(int i=0; i < _ni; ++i) {
//        for(int j=0; j < _nj; ++j) {
//            for(int k=0; k < _nk; ++k) {
//                _A[i][j][k] -= d;
//            }
//        }
//    }
//    return *this;
//}
Matrix& Matrix::operator*=(const Matrix &m)
{
    //for(int i=0; i < _ni; ++i) {
    //    for(int j=0; j < _nj; ++j) {
    //        for(int k=0; k < _nk; ++k) {
    //            _A[i][j][k] *= m._A[i][j][k];
    //        }
    //    }
    //}
    for(int idx=0; idx < _n; ++idx) {
        _A[idx] *= m._A[idx];
    }
    return *this;
}
//Matrix& Matrix::operator*=(const double d)
//{
//    for(int i=0; i < _ni; ++i) {
//        for(int j=0; j < _nj; ++j) {
//            for(int k=0; k < _nk; ++k) {
//                _A[i][j][k] *= d;
//            }
//        }
//    }
//    return *this;
//}
//
//
//Matrix operator-(const Matrix &a1, const Matrix &a2)
//{
//#ifdef DEBUG
//    assert(a1._ni==a2._ni);
//    assert(a1._nj==a2._nj);
//    assert(a1._nk==a2._nk);
//#endif
//    int ni = a1._ni;
//    int nj = a1._nj;
//    int nk = a1._nk;
//    Matrix A = Matrix(ni,nj,nk);
//    for(int i=0; i < ni; ++i) {
//        for(int j=0; j < nj; ++j) {
//            for(int k=0; k < nk; ++k) {
//                A(i,j,k) = a1._A[i][j][k] - a2._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
//Matrix operator-(const Matrix &m, const double d)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = m._A[i][j][k] - d;
//            }
//        }
//    }
//    return A;
//}
//Matrix operator-(const double d, const Matrix &m)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = d - m._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
//
//Matrix operator*(const Matrix &a1, const Matrix &a2)
//{
//#ifdef DEBUG
//    assert(a1._ni==a2._ni);
//    assert(a1._nj==a2._nj);
//    assert(a1._nk==a2._nk);
//#endif
//    int ni = a1._ni;
//    int nj = a1._nj;
//    int nk = a1._nk;
//    Matrix A = Matrix(ni,nj,nk);
//    for(int i=0; i < ni; ++i) {
//        for(int j=0; j < nj; ++j) {
//            for(int k=0; k < nk; ++k) {
//                A(i,j,k) = a1._A[i][j][k] * a2._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
//Matrix operator*(const Matrix &m, const double d)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = m._A[i][j][k] * d;
//            }
//        }
//    }
//    return A;
//}
//Matrix operator*(const double d, const Matrix &m)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = d * m._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
//
//Matrix operator/(const Matrix &a1, const Matrix &a2)
//{
//#ifdef DEBUG
//    assert(a1._ni==a2._ni);
//    assert(a1._nj==a2._nj);
//    assert(a1._nk==a2._nk);
//#endif
//    int ni = a1._ni;
//    int nj = a1._nj;
//    int nk = a1._nk;
//    Matrix A = Matrix(ni,nj,nk);
//    for(int i=0; i < ni; ++i) {
//        for(int j=0; j < nj; ++j) {
//            for(int k=0; k < nk; ++k) {
//                A(i,j,k) = a1._A[i][j][k] / a2._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
//Matrix operator/(const Matrix &m, const double d)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = m._A[i][j][k] / d;
//            }
//        }
//    }
//    return A;
//}
//Matrix operator/(const double d, const Matrix &m)
//{
//    Matrix A = Matrix(m._ni,m._nj,m._nk);
//    for(int i=0; i < m._ni; ++i) {
//        for(int j=0; j < m._nj; ++j) {
//            for(int k=0; k < m._nk; ++k) {
//                A(i,j,k) = d / m._A[i][j][k];
//            }
//        }
//    }
//    return A;
//}
