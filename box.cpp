#include "main.h"
#include "box.h"
//#include "matrix.h"

Box2d::Box2d(int nx, int ny, double ds) 
    : _nx(nx), _ny(ny), _ds(ds)
{
    _p = new Matrix(nx,ny);
    _u = new Matrix(nx,ny);
    _v = new Matrix(nx,ny);
    if(verbose) std::cout << "Initialized " 
        << nx << " by " << ny 
        << " box grid with ds=" << ds << "\n";
}
Box2d::~Box2d()
{
    delete _p;
    delete _u;
    delete _v;
}

void Box2d::init(double p0, double u0, double v0)
{
    *_p = p0;
    *_u = u0;
    *_v = v0;
    if(verbose) std::cout << "set initial values: "
        << p0 << " " << u0 << " " << v0 << "\n";
}
