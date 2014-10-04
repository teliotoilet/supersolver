#include "main.h"
#include "box.h"

int main()
{
    if(verbose) std::cout 
        << "\nI see a bad-ass mother who don't take no crap off of nobody!\n\n";

    Box2d g = Box2d(51,21,0.1);
    g.init(0,1,0);

    std::cout << "m = m + m\n";
    g.u() = g.u() + g.u();
    //std::cout << "u: " << g.u() << "\n";

    std::cout << "m -= d\n";
    g.u() -= 0.25;
    //std::cout << "u: " << g.u() << "\n";

    std::cout << "m += m\n";
    g.u() += g.u();
    //std::cout << "u: " << g.u() << "\n";

    std::cout << "m *= d\n";
    g.u() *= 0.5;
    //std::cout << "u: " << g.u() << "\n";

    std::cout << "m += 2*m\n";
    g.u() += 2*g.u();
    //std::cout << "u: " << g.u() << "\n";
 
    std::cout << "m += m + m\n";
    g.u() += g.u() + g.u();
    //std::cout << "u: " << g.u() << "\n";
}
