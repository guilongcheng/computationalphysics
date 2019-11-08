#ifndef DIFFERENTIATION
#define DIFFERENTIATION
#include "math.h"
#include "stdio.h"
#include <iostream>
using namespace std;

namespace GLC
{
template <typename T>
class Differentiation
{
public:
    T *x_, *y_, *dy_, *ddy_;
    T h_;
    T(*pf_)
    (T);
    int n_;

    Differentiation(T *x, int n)
    {
        n_ = n;
        x_ = x;
        h_ = *(x_ + 1) - *x_;
        dy_ = new T[n]();
        ddy_ = new T[n]();
    }

    ~Differentiation()
    {
        delete[] dy_;
        delete[] ddy_;
    }

    void ThreePoint(T *y)
    {
        y_ = y;
        for (int i = 1; i < n_ - 1; i++)
        {
            *(dy_ + i) = ((*(y_ + i + 1)) - (*(y_ + i - 1))) / (2 * h_);
            *(ddy_ + i) = ((*(y_ + i + 1)) - 2 * (*(y_ + i)) + (*(y_ + i - 1))) / (h_ * h_);
        }
    }

    void ThreePoint(T (*pf)(T))
    {
        pf_ = pf;
        for (int i = 1; i < n_ - 1; i++)
        {
            *(dy_ + i) = (pf_(*(x_ + i + 1)) - pf_(*(x_ + i - 1))) / (2 * h_);
            *(ddy_ + i) = (pf_(*(x_ + i + 1)) - 2 * pf_(*(x_ + i)) + pf_(*(x_ + i - 1))) / (h_ * h_);
        }
    }

    void adaptive(T (*pf)(T), double del = 1e-6, int maxstep = 1000)
    {
        int step;
        T d_, d1_, h1_;
        pf_ = pf;
        for (int i = 0; i < n_; i++)
        {
            d_ = (pf_(*(x_ + i) + h_) - pf_(*(x_ + i) - h_)) / (2 * h_);
            h1_ = h_ / 2;
            d1_ = (pf_(*(x_ + i) + h1_) - pf_(*(x_ + i) - h1_)) / (2 * h1_);
            step = 1;
            while ((h1_ * h1_ * abs(d_ - d1_) > del) && (step < maxstep))
            {
                h1_ = h1_ / 2;
                d_ = d1_;
                d1_ = (pf_(*(x_ + i) + h1_) - pf_(*(x_ + i) - h1_)) / (2 * h1_);
                step = step + 1;
            }
            *(dy_ + i) = (4 * d1_ - d_) / 3;
            cout << "The point " << i << " dy is " << *(dy_ + i) << ", Total iter is " << step << endl;
        }
    }

    void show()
    {
        for (int i = 1; i < n_ - 1; i++)
        {
            cout << i << ", dy = " << *(dy_ + i) << ", ddy = " << *(ddy_ + i) << endl;
        }
    }
};
} // namespace GLC
#endif