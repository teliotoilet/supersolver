#ifndef BOX_H
#define BOX_H

#include "matrix.h"

class Box2d
{
private:
    int _nx, _ny;
    double _ds;
    Matrix* _p;
    Matrix* _u;
    Matrix* _v;
    Matrix* _F;
    Matrix* _G;

    double _pinf, _uinf, _vinf;

    void calcXFlux();
    void calcYFlux();
    void calcXJac();
    void calcYJac();

public:
    Box2d(int nx, int ny, double ds=0.1); // nx,ny are cell counts
    ~Box2d();

    int Nx() const { return _nx; }
    int Ny() const { return _ny; }
    double ds() const { return _ds; }

    Matrix& p() { return *_p; }
    Matrix& u() { return *_u; }
    Matrix& v() { return *_v; }

    void setFreestream(double p0, double u0, double v0);
    void init(double p0, double u0, double v0);
    void initFS();

    void updateConvectiveFlux();
    void enforceConvectiveBCs();

};

#endif
