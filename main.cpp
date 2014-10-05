#include "main.h"
#include "box.h"

int main()
{
    if(verbose) std::cout 
        << "\nI see a bad-ass mother who don't take no crap off of nobody!\n\n";

    Box2d g = Box2d(51,21,0.1);
    g.setFreestream(0.0,1.0,0.0);
    g.initFS();

    // main loop
    g.updateConvectiveFlux();

}
