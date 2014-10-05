#ifndef MATRIX_H
#define MATRIX_H

#include <iostream>
#include <vector>

class Matrix
{
private:
    std::vector<double> _A;
    size_t _ni, _nj, _nk, _n;

public:
    Matrix(int ni, int nj=1, int nk=1);
    ~Matrix();

    int Ni() const { return _ni; }
    int Nj() const { return _nj; }
    int Nk() const { return _nk; }

    double sum();

    friend std::ostream& operator<< (std::ostream &out, Matrix &m);

    double& operator() (const int i);
    double& operator() (const int i, const int j);
    double& operator() (const int i, const int j, const int k);

    Matrix& operator= (const Matrix &m);
    Matrix& operator= (const double d);
    Matrix& operator= (const double d[]);
    //Matrix& operator+= (const Matrix &m);
    //Matrix& operator+= (const double d);
    //Matrix& operator-= (const Matrix &m);
    //Matrix& operator-= (const double d);
    Matrix& operator*= (const Matrix &m);
    //Matrix& operator*= (const double d);

    //friend Matrix operator+ (const Matrix &a1, const Matrix &a2);
    //friend Matrix operator+ (const Matrix &m, const double d);
    //friend Matrix operator+ (const double d, const Matrix &m);

    //friend Matrix operator- (const Matrix &a1, const Matrix &a2);
    //friend Matrix operator- (const Matrix &m, const double d);
    //friend Matrix operator- (const double d, const Matrix &m);

    //friend Matrix operator* (const Matrix &a1, const Matrix &a2);
    //friend Matrix operator* (const Matrix &m, const double d);
    //friend Matrix operator* (const double d, const Matrix &m);

    //friend Matrix operator/ (const Matrix &a1, const Matrix &a2);
    //friend Matrix operator/ (const Matrix &m, const double d);
    //friend Matrix operator/ (const double d, const Matrix &m);

};

#endif
