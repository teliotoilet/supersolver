#include "main.h"
#include "matrix.h"

#include <cassert>

Matrix::Matrix(int ni, int nj, int nk)
    : _ni(ni), _nj(nj), _nk(nk)
{
    _A = new double**[_ni];
    for(int i=0; i < _ni; ++i) {
        _A[i] = new double*[_nj];
        for(int j=0; j < _nj; ++j) {
            _A[i][j] = new double[_nk];
        }
    }
    if(verbose) std::cout << "Constructed matrix with dimensions (" 
        << ni << "," << nj << "," << nk << ")\n";
}
//Matrix::Matrix(const Matrix &m)
//    : _ni(m._ni), _nj(m._nj), _nk(m._nk)
//{
//    for(int i=0; i < _ni; ++i) {
//        for(int j=0; j < _nj; ++j) {
//            for(int k=0; k < _nk; ++k) {
//                _A[i][j][k] = m._A[i][j][k];
//            }
//        }
//    }
//    if(verbose) std::cout << "Copied matrix with dimensions (" 
//        << m._ni << "," << m._nj << "," << m._nk << ")\n";
//}
Matrix::~Matrix()
{
    if(_ni > 0) {
        for(int i=0; i < _ni; ++i) {
            for(int j=0; j < _nj; ++j) {
                delete[] _A[i][j];
            }
            delete[] _A[i];
        }
        delete[] _A;
    }
    if(verbose) std::cout << "destroyed matrix with dimensions ("
        << _ni << "," << _nj << "," << _nk << ")\n";
}

std::ostream& operator<< (std::ostream &out, Matrix &m)
{
    out << "[";
    for(int k=0; k < m._nk; ++k) {

        if(k==0) out << "[";
        else out << ",";
        for(int j=0; j < m._nj; ++j) {

            if(j==0) out << "[";
            else out << ",";
            for(int i=0; i < m._ni; ++i) {

                if(i==0) out << "[";
                else out << ",";

                out << m._A[i][j][k];

            }
            out << "]";
        }
        out << "]";
    }
    out << "]";
}

Matrix& Matrix::operator= (const Matrix &m)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] = m._A[i][j][k];
            }
        }
    }
    if(verbose) std::cout << "Overloaded assignment with dimensions (" 
        << m._ni << "," << m._nj << "," << m._nk << ")\n";
    return *this;
}
Matrix& Matrix::operator= (const double d)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] = d;
            }
        }
    }
    return *this;
}

double& Matrix::operator() (const int i, const int j, const int k)
{
    assert(i >= 0 && i < _ni);
    assert(j >= 0 && j < _nj);
    assert(k >= 0 && k < _nk);
    return _A[i][j][k];
}

Matrix operator+(const Matrix &a1, const Matrix &a2)
{
    int ni = a1._ni; assert(a1._ni==a2._ni);
    int nj = a1._nj; assert(a1._nj==a2._nj);
    int nk = a1._nk; assert(a1._nk==a2._nk);
    Matrix A = Matrix(ni,nj,nk);
    for(int i=0; i < ni; ++i) {
        for(int j=0; j < nj; ++j) {
            for(int k=0; k < nk; ++k) {
                A(i,j,k) = a1._A[i][j][k] + a2._A[i][j][k];
            }
        }
    }
    return A;
}
Matrix operator+(const Matrix &m, const double d)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = m._A[i][j][k] + d;
            }
        }
    }
    return A;
}
Matrix operator+(const double d, const Matrix &m)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = d + m._A[i][j][k];
            }
        }
    }
    return A;
}
Matrix& Matrix::operator+=(const Matrix &m)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] += m._A[i][j][k];
            }
        }
    }
    return *this;
}
Matrix& Matrix::operator+=(const double d)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] += d;
            }
        }
    }
    return *this;
}
Matrix& Matrix::operator-=(const Matrix &m)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] -= m._A[i][j][k];
            }
        }
    }
    return *this;
}
Matrix& Matrix::operator-=(const double d)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] -= d;
            }
        }
    }
    return *this;
}
Matrix& Matrix::operator*=(const Matrix &m)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] *= m._A[i][j][k];
            }
        }
    }
    return *this;
}
Matrix& Matrix::operator*=(const double d)
{
    for(int i=0; i < _ni; ++i) {
        for(int j=0; j < _nj; ++j) {
            for(int k=0; k < _nk; ++k) {
                _A[i][j][k] *= d;
            }
        }
    }
    return *this;
}


Matrix operator-(const Matrix &a1, const Matrix &a2)
{
    int ni = a1._ni; assert(a1._ni==a2._ni);
    int nj = a1._nj; assert(a1._nj==a2._nj);
    int nk = a1._nk; assert(a1._nk==a2._nk);
    Matrix A = Matrix(ni,nj,nk);
    for(int i=0; i < ni; ++i) {
        for(int j=0; j < nj; ++j) {
            for(int k=0; k < nk; ++k) {
                A(i,j,k) = a1._A[i][j][k] - a2._A[i][j][k];
            }
        }
    }
    return A;
}
Matrix operator-(const Matrix &m, const double d)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = m._A[i][j][k] - d;
            }
        }
    }
    return A;
}
Matrix operator-(const double d, const Matrix &m)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = d - m._A[i][j][k];
            }
        }
    }
    return A;
}

Matrix operator*(const Matrix &a1, const Matrix &a2)
{
    int ni = a1._ni; assert(a1._ni==a2._ni);
    int nj = a1._nj; assert(a1._nj==a2._nj);
    int nk = a1._nk; assert(a1._nk==a2._nk);
    Matrix A = Matrix(ni,nj,nk);
    for(int i=0; i < ni; ++i) {
        for(int j=0; j < nj; ++j) {
            for(int k=0; k < nk; ++k) {
                A(i,j,k) = a1._A[i][j][k] * a2._A[i][j][k];
            }
        }
    }
    return A;
}
Matrix operator*(const Matrix &m, const double d)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = m._A[i][j][k] * d;
            }
        }
    }
    return A;
}
Matrix operator*(const double d, const Matrix &m)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = d * m._A[i][j][k];
            }
        }
    }
    return A;
}

Matrix operator/(const Matrix &a1, const Matrix &a2)
{
    int ni = a1._ni; assert(a1._ni==a2._ni);
    int nj = a1._nj; assert(a1._nj==a2._nj);
    int nk = a1._nk; assert(a1._nk==a2._nk);
    Matrix A = Matrix(ni,nj,nk);
    for(int i=0; i < ni; ++i) {
        for(int j=0; j < nj; ++j) {
            for(int k=0; k < nk; ++k) {
                A(i,j,k) = a1._A[i][j][k] / a2._A[i][j][k];
            }
        }
    }
    return A;
}
Matrix operator/(const Matrix &m, const double d)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = m._A[i][j][k] / d;
            }
        }
    }
    return A;
}
Matrix operator/(const double d, const Matrix &m)
{
    Matrix A = Matrix(m._ni,m._nj,m._nk);
    for(int i=0; i < m._ni; ++i) {
        for(int j=0; j < m._nj; ++j) {
            for(int k=0; k < m._nk; ++k) {
                A(i,j,k) = d / m._A[i][j][k];
            }
        }
    }
    return A;
}
