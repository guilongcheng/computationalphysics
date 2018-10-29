#include <iostream>
#include "../Eigen/Core"
using namespace Eigen;
void print_size(const Ref<const MatrixXf>& b)
{
  std::cout << "size (rows, cols): " << b.size() << " (" << b.rows()
            << ", " << b.cols() << ")" << std::endl;
  std::cout<<b(0)<<std::endl;
}
int main()
{
    Vector3f v;
//    print_size(v);
    // v.asDiagonal() returns a 3x3 diagonal matrix pseudo-expression
    Matrix3f M;
    M = v.asDiagonal();
    print_size(M);
}
