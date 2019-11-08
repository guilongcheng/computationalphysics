#!/bin/env python
import numpy as np

def Roberg(f,a,b,eps=10**(-6),nmax=1000,flag="",**args):
    """龙贝格积分

    :f: TODO
    :a: TODO
    :b: TODO
    :eps: TODO
    :returns: TODO

    """
    M=1;tol=10
    k = 0
    T = np.zeros((nmax,nmax))
    h = b - a
    T[0,0] = (h/2) * (f(a,**args) + f(b,**args))
    while tol > eps:
        k = k + 1
        h = h / 2
        x = a + h * (2 * np.arange(1,M+1) - 1)
        Q = np.sum(f(x,**args))
        T[k,0] = T[k-1,0] / 2 + h*Q
        M = 2*M
        for j in range(k):
            T[k,j+1] = T[k,j] + (T[k,j] - T[k-1,j]) / (4**(j+1) - 1)
        tol = abs(T[k,k] - T[k-1,k-1])
    print("step is ",k,"  龙贝格积分结果为:",T[k,k])
    if "show" in flag:        
        print("龙贝格表",T[:k+1,:k+1])
    return T[k,k]

def main():
    print("==========例题4.2.7==========")
    def f(x):
        return np.log(1+x)/(1+x**2)
    T = Roberg(f,0,1,10**(-6))
    print("龙贝格积分结果为:",T)

    print("==========例题4.2.8==========")
    def f(x):
        return 4/(1+x**2)
    T = Roberg(f,0,1,10**(-6))
    print("龙贝格积分结果为:",T)

if __name__ == "__main__":
    main()



