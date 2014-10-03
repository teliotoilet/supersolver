#include "main.h"
#include "box.h"

int main()
{
    if(verbose) std::cout 
        << "\nI see a bad-ass mother who don't take no crap off of nobody!\n\n";

    Box2d g = Box2d(51,21,0.1);
    std::cout << "u: " << g.u() << "\n";
    g.init(0,1,0);
    std::cout << "u: " << g.u() << "\n";

    std::cout << "manipulating u\n";
    g.u() = g.u() + g.u();
    std::cout << "u: " << g.u() << "\n";

    std::cout << "manipulating u second time\n";
    g.u() = 0.5 + g.u();
    g.u() = g.u() - 0.25;
    std::cout << "u: " << g.u() << "\n";
}
