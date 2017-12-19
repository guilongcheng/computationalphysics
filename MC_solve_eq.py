#!/bin/env python

import numpy as np

def MCSolveEq1(f,a,b,N):
    xi = np.random.random(N)
    eta = np.zeros(len(xi))
    for i in range(len(xi)):
        if f(a+(b-a)*xi[i]) >= 0:
            eta[i]=1
    if f(a)<0:
        root = a+(b-a)*np.average(eta)
    if f(a)>0:
        root = a+(b-a)*(1-np.average(eta))
    print("root is ",root)

def MCSolveEq2(f,a,b,N):
    xi = np.random.random(N)
    imax=0.0
    imin=0.0
    ii = 0
    ij = 0
    for i in np.nditer(xi):
        if f(a)<0:
            if f(a-(b-a)*i) <= 0  :
                ii=ii+1
                if ii == 1:
                    imax = i
                if f(a-(b-a)*i) > f(a-(b-a)*imax):
                    imax = i
            if f(a-(b-a)*i) >= 0 and f(a-(b-a)*i) < f(a-(b-a)*imin):
                ij=ij+1
                if ij == 1:
                    imin = i
                if f(a-(b-a)*i) > f(a-(b-a)*imin):
                    imin = i
        if f(a) > 0:
            if f(a-(b-a)*i) >= 0  :
                ii=ii+1
                if ii == 1:
                    imax = i
                if f(a-(b-a)*i) > f(a-(b-a)*imax):
                    imax = i
            if f(a-(b-a)*i) <= 0 and f(a-(b-a)*i) < f(a-(b-a)*imin):
                ij=ij+1
                if ij == 1:
                    imin = i
                if f(a-(b-a)*i) > f(a-(b-a)*imin):
                    imin = i
    root = (imax + imin) / 2.0
    print("root is ", root )
    return root

def main():
    def f(x):
        return x**4 -10*x**3 +35*x**2 -50*x +24
    def f2(x):
        return x**2 - 2*x + 1

    find1 = MCSolveEq1(f,0.7,1.2,10000)
    find2 = MCSolveEq2(f,0.7,1.2,10000)
    find1 = MCSolveEq1(f2,0.7,1.2,10000)
    find2 = MCSolveEq2(f2,0.7,1.2,10000)

if __name__ == "__main__":
    main()


