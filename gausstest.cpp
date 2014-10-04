#include "main.h"
#include "matrix.h"
#include "util.h"

int main()
{
    double lhs[] = { 10, -1,  2,  0,
                     -1, 11, -1,  3,
                      2, -1, 10, -1,
                      0,  3, -1,  8};
    double rhs[] = { 6, 25, -11, 15};
    double guess[] = { 0, 0, 0, 0 };
    Matrix A = Matrix(4,4,lhs);
    Matrix b = Matrix(4,rhs);
    Matrix x = Matrix(4,guess);


    // solver iterations
    double R,Rlast = 0.0;
    for(int n=0; n<10; ++n) {

        for(int i=0; i < 4; ++i) {
            double s = 0.0;
            for(int j=0; j < 4; ++j) {
                if(j != i) {
                    s += A(i,j)*x(j);
                }
            }
            x(i) = (b(i)-s)/A(i,i);
        }
        R = calcResid(x,A,b);
        std::cout << "iter " << n
            << ": x= " 
            << x(0) << " "
            << x(1) << " "
            << x(2) << " "
            << x(3) << " "
            << "residuals: "
            << R << " " << R/Rlast
            << std::endl;

        Rlast = R;
    }

}
