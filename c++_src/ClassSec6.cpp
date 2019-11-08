#include "SolveRoot.h"
#include "Differentiation.h"
#include <iostream>
#include "math.h"
#define PI 3.14159265358979323846264338327950288419716939937510
using namespace std;
using namespace GLC;

double func(double x)
{
    return exp(x) * log(x) - x * x;
}

double dfunc(double x)
{
    return exp(x) *(log(x) + 1.0 / x) - 2*x;
}

int main(int argc, char const *argv[])
{
    cout << "计算微分" << endl;
    cout << "====================================================" << endl;
    double *x, *y;
    double h;
    int n = 6;
    x = new double[n];
    y = new double[n];

    for (int i = 0; i < n; i++)
    {
        *(x + i) = i * PI / (2 * (n - 1));
        *(y + i) = sin((*(x + i)));
        // cout<<"x[i] = "<<*(x + i)<<"; y[i] = "<<*(y+i)<<endl;
    }
    h = (*(x + 1)) - (*(x));
    cout << "h is " << h << endl;
    cout << "Three point results" << endl;
    Differentiation<double> diff(x, n);
    diff.ThreePoint(y);
    diff.show();

    for (int i = 1; i < n - 1; i++)
    {
        cout << " cos(x) = " << cos(*(x + i)) << endl;
    }

    cout << "adaptive scheme results" << endl;
    diff.adaptive(cos);

    cout << "方程求根" << endl;
    cout << "====================================================" << endl;
    SolveRoot<double> Solve(func);
    Solve.Bisection(1,2);
    Solve.Newton(1,dfunc);
    Solve.Secant(1);
    return 0;
}
