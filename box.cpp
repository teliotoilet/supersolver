#include "main.h"
#include "box.h"
//#include "matrix.h"

Box2d::Box2d(int nx, int ny, double ds)
    : _nx(nx), _ny(ny), _ds(ds)
{
    _p = new Matrix(nx,ny);     // cell centered
    _u = new Matrix(nx+1,ny+1); // node centered
    _v = new Matrix(nx+1,ny+1); // node centered
    _F = new Matrix(nx,ny,2);   // interior edges
    _G = new Matrix(nx,ny,2);   // interior edges
    if(verbose) std::cout << "Initialized box grid with "
        << nx << " by " << ny 
        << " cells, ds=" << ds << "\n";
}
Box2d::~Box2d()
{
    delete _p;
    delete _u;
    delete _v;
}

void Box2d::setFreestream(double p0, double u0, double v0)
{
    _pinf = p0;
    _uinf = u0;
    _vinf = v0;
}

void Box2d::initFS()
{
    init(_pinf,_uinf,_vinf);
}

void Box2d::init(double p0, double u0, double v0)
{
    *_p = p0;
    *_u = u0;
    *_v = v0;
    if(verbose) std::cout << "set initial values: "
        << p0 << " " << u0 << " " << v0 << "\n";
}

void Box2d::updateConvectiveFlux()
{
    calcXFlux();
    calcYFlux();
}

void Box2d::enforceConvectiveBCs()
{
}

void Box2d::calcXFlux()
{
    double pi, ui, vi;

    // this is a really roundabout way of doing this...
    Matrix p = *_p;
    Matrix u = *_u;
    Matrix v = *_v;
    Matrix F = *_F;

    for(int i=0; i < _nx; ++i) {
        for(int j=0; j < _ny; ++j) {

            // average p from cell centers to edge
            if(j==0) pi = _pinf;
            else pi = (p(i,j-1) + p(i,j))/2.0;

            // averge velocities from nodes to edge
            ui = (u(i,j) + u(i+1,j))/2.0;
            vi = (v(i,j) + v(i+1,j))/2.0;

            // calculate flux
            F(i,j,0) = ui*ui + pi;
            F(i,j,1) = ui*vi;

        }
    }

    *_F = F;
}

void Box2d::calcYFlux()
{
    double pi, ui, vi;

    // this is a really roundabout way of doing this...
    Matrix p = *_p;
    Matrix u = *_u;
    Matrix v = *_v;
    Matrix G = *_F;

    for(int i=0; i < _nx; ++i) {
        for(int j=0; j < _ny; ++j) {

            // average p from cell centers to edge
            if(i==0) pi = _pinf;
            else pi = (p(i-1,j) + p(i,j))/2.0;

            // averge velocities from nodes to edge
            ui = (u(i,j) + u(i,j+1))/2.0;
            vi = (v(i,j) + v(i,j+1))/2.0;

            // calculate flux
            G(i,j,1) = ui*vi;
            G(i,j,2) = vi*vi + pi;

        }
    }

    *_G = G;
}
