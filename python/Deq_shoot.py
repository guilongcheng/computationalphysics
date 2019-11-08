#!/bin/env python
import numpy as np
import scipy.integrate as spint
import matplotlib.pyplot as plt
from Deq_RK import RungKutta

def shoot(fun,x,y0,g1,eps=0.00001,**args):
    """TODO: Docstring for shoot.

    :fun: TODO
    :x: TODO
    :y0: TODO
    :g1: TODO
    :returns: TODO

    """
    test_y0=np.array([y0[0],g1])
    y = RungKutta(fun,x,test_y0,**args)
    yb=y[0,-1]

    test_y0[1] = y0[1]*g1/y[0,-1]
    g2 = test_y0[1]
    y = RungKutta(fun,x,test_y0,**args)

    delta = y[0,-1] - y0[1]
    while abs(delta) > eps:
        print(g1,g2,(y0[1]-y[0,-1]),g2-g1,y[0,-1]-yb)
        test_y0[1] = g2+(y0[1]-y[0,-1])* \
            (g2-g1)/(y[0,-1]-yb)
        g1 = g2
        g2 = test_y0[1]
        yb = y[0,-1]
        y = RungKutta(fun,x,test_y0,**args)
        delta = y[0,-1] - y0[1]

    return y


def main():
    print("==========例题5.3.3==========")
    def fun1(x,y):
        return np.array([y[1],-9.8-0.01*y[1]])
    def fun2(x,y):
        return np.array([y[1],-9.8])

    x = np.linspace(0,5,100)
    y0=np.array([0,40])
    g1 = 30.0 # 注意需要给定为实数。即加上.0的后缀。

    y1 = shoot(fun1,x,y0,g1)
    print("v0 = ",y1[1,0])

    y2 = shoot(fun2,x,y0,g1)
    print("v0 = ",y2[1,0])

    def bc(ya,yb):
        return np.array([ya[0],yb[0]-40]) #注意边界条件的写法。
    p0 = np.zeros((2,x.size)) #这里的p0并不重要，初始值体现在边界条件里。
    y3 = spint.solve_bvp(fun1,bc,x,p0)
    print("v0 = ",y3.sol(x)[1,0])

    plt.xlabel("t")
    plt.ylabel("y")
    plt.plot(x,y1[0],'-',label="f on")
    plt.plot(x,y2[0],'-.',label="f off")
    plt.plot(x,y3.sol(x)[0],'o',label='sp')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()


