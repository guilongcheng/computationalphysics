#ifndef SOLVEROOT
#define SOLVEROOT
#include <iostream>
#include <math.h>
using namespace std;
namespace GLC
{
template <typename T>
class SolveRoot
{
private:
    /* data */
public:
    T(*pf_)
    (T);
    T solve;
    SolveRoot(T (*pf)(T))
    {
        pf_ = pf;
    };

    ~SolveRoot(){};

    void Bisection(T a, T b, double del = 1e-6, int maxstep = 1000)
    {
        T f1, f2;
        T x1, x2, x;
        int step;
        x1 = a;
        x2 = b;
        step = 0;
        while (abs(x1 - x2) > del && step < maxstep)
        {
            x = (x1 + x2) / 2;
            f1 = pf_(x1);
            f2 = pf_(x);
            if (f1 * f2 < 0)
            {
                x2 = x;
            }
            else
            {
                x1 = x;
            }
            step++;
        }
        if (step >= maxstep)
        {
            cout << "Get the max iters " << step << endl;
        }
        solve = (x1 + x2) / 2.0;
        cout << "Total steps are " << step << "; The bisection solve is " << solve << endl;
    }

    void Newton(T x0, T (*pdf)(T), double del = 1e-6, int maxstep = 1000)
    {
        int step = 0;
        T f, df;
        T x, dx;
        x = x0;
        f = pf_(x);
        df = pdf(x);
        dx = -f / df;
        while (abs(dx) > del && step < maxstep)
        {
            x = x + dx;
            f = pf_(x);
            df = pdf(x);
            dx = -f / df;
            step++;
        }
        if (step >= maxstep)
        {
            cout << "Get the max iters " << step << endl;
        }
        solve = x + dx;
        cout << "Total steps are " << step << "; The Newton solve is " << solve << endl;
    }

    void Secant(T x0,T x01 = 0.1, double del = 1e-6, int maxstep = 1000)
    {
        int step = 0;
        T f, f1 ;
        T x1, x2, x, dx;
        x1 = x0;
        x = x0 + x01;
        f1 = pf_(x1);
        f = pf_(x);
        dx = -f * (x - x1) / (f - f1);
        while (abs(x1 - x) > del && step < maxstep)
        {
            x1 = x;
            x = x + dx;
            f1 = pf_(x1);
            f = pf_(x);
            dx = -f * (x - x1) / (f - f1);
            step++;
        }
        if (step >= maxstep)
        {
            cout << "Get the max iters " << step << endl;
        }
        solve = x + dx;
        cout << "Total steps are " << step << "; The Secant solve is " << solve << endl;
    }
};
} // namespace GLC
#endif