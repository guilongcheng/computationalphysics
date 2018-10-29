#include"iostream"
#include"../Eigen/Dense"
#include"lib_linalg.h"
using namespace std;
using namespace Eigen;
int main()
{
    Matrix3d m;
    Vector3d b;
    Vector3d x;
    m<<1,2,1,
       2,2,3,
       -1,0,-3;
    b<<0,3,0;

    gauss_jordan(m,b,x);

    showM(m);
    showV(x);
    return 0;
}
