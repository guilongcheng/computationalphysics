#ifndef DSOLVE
#define DSOLVE
#include <iostream>
#include <math.h>
using namespace std;
namespace GLC
{
template <typename T>
class DSolve
{
private:
    /* data */
public:
    T(*pf)
    (T,T);
    T *y;
    int n;
    T solve;
    DSolve(T y0_, T (*pf_)(T,T), int n_ )
    {
        pf = pf_;
        n = n_;
        y = new T[n]();
        *y = y0_;
    };

    ~DSolve(){
        delete[] y;
    };

    void Eular(T t0, T tau){
        for (int i = 0; i < n-1; i++)
        {
           *(y + i + 1) = *(y + i) + tau * pf(*(y + i),t0 + i * tau ); 
        }
    }

    void PreEular(T t0, T tau){
        for (int i = 0; i < n-1; i++)
        {
           *(y + i + 1) = *(y + i) + tau * pf(*(y + i),t0 + i * tau ); 
           *(y + i + 1) = *(y + i) + tau / 2.0 * (pf(*(y + i),t0 + i * tau ) + pf(*(y + i + 1),t0 + (i + 1) * tau )); 
        }
    }

    void RungeKutta1(T t0, T tau){
        T c1,c2,c3,c4;
        for (int i = 0; i < n-1; i++)
        {
            c1 = tau * pf(*(y + i), t0 + i * tau);
            c2 = tau * pf(*(y + i) + c1 / 2, t0 + i * tau + tau / 2);
            c3 = tau * pf(*(y + i) + c2 / 2, t0 + i * tau + tau / 2);
            c4 = tau * pf(*(y + i) + c3 , t0 + i * tau + tau );
            *(y + i + 1) = *(y + i) + 1.0 / 6.0 * (c1 + 2*c2 + 2*c3 + c4) ;
        }
    }

    void RungeKutta2(T t0, T tau){
        T c1,c2,c3,c4;
        for (int i = 0; i < n-1; i++)
        {
            c1 = tau * pf(*(y + i), t0 + i * tau);
            c2 = tau * pf(*(y + i) + c1 / 2, t0 + i * tau + tau / 2);
            c3 = tau * pf(*(y + i) + c2 / 2, t0 + i * tau + tau / 2);
            c4 = tau * pf(*(y + i) + c3 , t0 + i * tau + tau );
            *(y + i + 1) = *(y + i) + 1.0 / 6.0 * (c1 + 2*c2 + 2*c3 + c4) ;
        }
    }

    void show(){
        for (int i = 0; i < n; i++)
        {
            cout<<"y("<<i<<") = "<<*(y + i)<<endl;
        }
        
    }

};
} // namespace GLC
#endif