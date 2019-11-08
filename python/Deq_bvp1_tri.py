#!/bin/env python
import numpy as np
import scipy.integrate as spint
from linalg_tridag import tridag
import matplotlib.pyplot as plt


def bvp1_tri(u,v,w,f,x,p0):
    h = x[1] - x[0]
    a = u - 0.5*h*v
    b = h**2*w -2*u
    c = u + 0.5*h*v
    d = h**2 * f
    a[0] = 0;b[0] = 1;c[0] = 0; d[0] = p0[0]
    a[-1] = 0;b[-1] =1;c[-1] =0;d[-1] = p0[1]
    y = tridag(a,b,c,d)
    return y

def bvp2_tri(u,v,w,f,x,p0):
    h = x[1] - x[0]
    a = u - 0.5*h*v
    b = h**2*w -2*u
    c = u + 0.5*h*v
    d = h**2 * f
    a[0] = 0;b[0] = 1;c[0] = 0; d[0] = p0[0]
    a[-1] = -1;b[-1] = 1;c[-1] =0;d[-1] = h*p0[1]
    y = tridag(a,b,c,d)
    return y

def main():
    print("==========例题5.3.1==========")
    #  利用三对角矩阵求解
    p0 = [0,0]
    x = np.linspace(2,3,11)
    n = len(x)
    u = -1*np.ones(n)
    v = np.zeros(n)
    w = 2/x**2
    f = 1/x
    y = bvp1_tri(u,v,w,f,x,p0)

    # 调用scipy.integrate 的函数solve_bvp求解
    def fun(x,y):
        return [y[1],-1/x+2/x**2*y[0]]
    def bc(ya,yb):
        return np.array([ya[0],yb[0]])
    y0 = np.zeros((2,n))
    y2 = spint.solve_bvp(fun,bc,x,y0)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(x,y,'-',label='tri')
    plt.plot(x,y2.sol(x)[0],'o',label="sp")
    plt.grid()
    plt.legend()
    plt.show()

    print("==========例题5.3.2==========")

    #  利用三对角矩阵求解
    p0 = [0,0]
    x = np.linspace(0,2,11)
    n = len(x)
    u = np.ones(n)
    v = np.zeros(n)
    w = 9*np.ones(n)
    f = x
    y = bvp2_tri(u,v,w,f,x,p0)

    # 调用scipy.integrate 的函数solve_bvp求解
    def fun(x,y):
        return [y[1],x-9*y[0]]
    def bc(ya,yb):
        return np.array([ya[0],yb[1]])# 这里yb[1] 代表着y'(b)
    y0 = np.zeros((2,n))
    y2 = spint.solve_bvp(fun,bc,x,y0)

    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(x,y,'-',label='tri')
    plt.plot(x,y2.sol(x)[0],'o',label='sp')
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()

