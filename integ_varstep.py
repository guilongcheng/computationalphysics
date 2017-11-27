#!/bin/env python
import numpy as np
import integ_Romberg

def NaivTrap(f,a,b,eps):
    """变步长梯形积分

    :f: TODO
    :a: TODO
    :b: TODO
    :eps: TODO
    :returns: TODO

    """
    n = 1
    h = (b-a)/n
    I1 = 0.5*h*(f(a)+f(b))
    tol=1;
    while tol>eps:
        I0 = I1
        n = 2*n
        h = (b-a)/n
        x = np.linspace(a,b,n+1)
        y = f(x)
        I1 = h*(0.5*y[0] + np.sum(y[1:n]) +  0.5*y[n])
        tol=abs(I1-I0)

    print("n is ",n)
    return I1

def aTrapInt(f,a,b,eps):
    """变步长梯形积分

    :f: TODO
    :a: TODO
    :b: TODO
    :eps: TODO
    :returns: TODO

    """
    tol=1;nsub=1
    inall = 0
    T = 0.5*(b-a)*(f(a)-f(b))
    while tol>eps:
        T0=T
        nsub = 2*nsub
        n = nsub+1
        h = (b-a) / nsub
        x = np.linspace(a,b,n)
        inall = inall + np.sum(f(x[1:n-1:2]))
        T = 0.5*h*(f(a) + 2*inall + f(b))
        tol = abs(T-T0)

    print("n is ",nsub)
    return T

def NaivSimp(f,a,b,eps):
    """变步长辛普森积分

    :f: TODO
    :a: TODO
    :b: TODO
    :eps: TODO
    :returns: TODO

    """
    n=2
    h = (b-a)/n
    I1 = h*(f(a) + 4*f(a+h) + f(b)) / 3
    tol = 1
    while tol > eps:
        I0 = I1
        n = 2*n
        h = (b-a)/n
        x = np.linspace(a,b,n+1)
        y = f(x)
        I1 = h * (y[0] + 4*np.sum(y[1:n+1:2]) + 2*np.sum(y[2:n:2]) + y[n])
        tol = abs(I1 - I0)

    return I1

def aSimpInt(f,a,b,eps):
    """变步长辛普森积分

    :f: TODO
    :a: TODO
    :b: TODO
    :eps: TODO
    :returns: TODO

    """
    nsub = 2
    odd = f((a+b)/2)
    even = 0
    S = (b-a) * (f(a) + 4*odd + 2*even + f(b)) / 6
    inall = even + odd
    tol = 1
    while tol > eps:
        S0 = S
        nsub = 2 * nsub
        n = nsub + 1
        h = (b - a) / nsub
        x = np.linspace(a,b,n)
        odd = np.sum(f(x[1:n:2]))
        even = inall
        S = (h/3) * (f(a) + 4*odd + 2*even + f(b))
        inall = even + odd
        tol = abs(S-S0)
    print("n is ",nsub)
    return S

def main():

    print("==========例题4.2.6==========")
    def f(x):
        return np.exp(-x**2)

    print("变步长梯形公式的结果为：",aTrapInt(f,0,1,10**(-6)))
    print("变步长辛普森公式的结果为：",aSimpInt(f,0,1,10**(-6)))

    print("==========作乐4.2.6==========")

    def f(x):
        return 4/(1+x**2)

    print("变步长梯形公式的结果为：",aTrapInt(f,0,1,10**(-5)))
    print("变步长辛普森公式的结果为：",aSimpInt(f,0,1,10**(-5)))
    T,k=integ_Romberg.Roberg(f,0,1,10**(-5))
    print("Table:")
    print(T[:k+1,:k+1])



if __name__ == "__main__":
    main()

