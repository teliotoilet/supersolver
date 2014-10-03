#ifndef BOX_H
#define BOX_H

#include "matrix.h"

class Box2d
{
private:
    const int _nx, _ny;
    const double _ds;
    Matrix* _p;
    Matrix* _u;
    Matrix* _v;

public:
    Box2d(int nx, int ny, double ds=0.1);
    ~Box2d();

    int Nx() const { return _nx; }
    int Ny() const { return _ny; }
    double ds() const { return _ds; }

    Matrix& p() { return *_p; }
    Matrix& u() { return *_u; }
    Matrix& v() { return *_v; }

    void init(double p0, double u0, double v0);

};

#endif
