#!/bin/env python
import numpy as np
import scipy as sp
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


def rk4n(f,n,x,h,y):
    """
    f:函数
    n:方程组的数量
    x:x点
    h:步长
    y:y的初始值
    """
    x0=x;y0=y
    dy1 =f(x0,y)
    y = y0 + 0.5*h*dy1
    dy2 = f(x0+0.5*h,y)
    y = y0 + 0.5*h*dy2
    dy3 = f(x0+0.5*h,y)
    y = y0 + h*dy3
    dy4 = f(x0+h,y)
    dy = (dy1 + 2*(dy2+dy3) + dy4) / 6
    y = y0 + h*dy
    return y

def rk4s(f,y0,tspan,h):
    """
    f:函数
    tspan:区间
    y0:初始值
    h:步长
    """
    t0=tspan[0];tf=tspan[1]
    t = np.arange(t0,tf+h,h)
    n = len(t)
    y = np.zeros((len(y0),n))
    y[:,0] = y0
    for i in range(1,n):
        s1 = f(t[i-1],y[:,i-1])
        s2 = f(t[i-1] + h/2, y[:,i-1] + h* s1/2)
        s3 = f(t[i-1]+h/2, y[:,i-1] + h* s2/2)
        s4 = f(t[i-1]+h, y[:,i-1] + h* s3)
        y[:,i] = y[:,i-1] + h*(s1+2*s2+2*s3+s4)/6
    return t,y

def main():
    print("==========例题5.2.4==========")

    def f(t,y):
        return np.array([y[3],y[4],y[5],y[4],-y[3],0])

    n = 6; t0=0; h=0.1
    y = np.array([1,1,0,1,0,1])
    y1 = np.zeros((500,6))
    t = np.arange(t0,500,h)

    for i in range(500):
        y = rk4n(f,n,t[i],h,y)
        y1[i,:] = y

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(y1[:,0],y1[:,1],y1[:,2])
    plt.show()

    print("==========例题5.2.5==========")

    def f(t,y):
        return np.array([y[1],2*y[1]-2*y[0]+np.exp(2*t)*np.sin(t)])

    tspan=[0,1]; h=0.02
    y0 = np.array([-0.4,-0.6])

    t,y = rk4s(f,y0,tspan,h)

    fig = plt.figure()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.plot(t,y[0],'-.',label='y(1)=y(x)')
    plt.plot(t,y[1],'-',label="y(2) = dy/dx")
    plt.legend()
    plt.show()

    print("==========例题5.2.7==========")

    def f(t,y):
        return np.array([998*y[0]+1998*y[1],-999*y[0]-1999*y[1]])

    tspan = [0,0.02]; h = 0.02/400
    y0 = np.array([1,0])

    t,y = rk4s(f,y0,tspan,h)

    fig = plt.figure()
    plt.xlabel("t")
    plt.ylabel("y")
    plt.plot(t,y[0],'-.',label='x(t)')
    plt.plot(t,y[1],'-',label="y(t)")
    plt.legend()
    plt.show()

    tspan = [0.02,6]; h = (6-0.02)/10000
    y0 = np.array([y[0,-1],y[1,-1]])

    t,y = rk4s(f,y0,tspan,h)

    plt.xlabel("t")
    plt.ylabel("y")
    plt.plot(t,y[0],'-.',label='x(t)')
    plt.plot(t,y[1],'-',label="y(t)")
    plt.legend()
    plt.show()

    tspan = [0.0,6]; h = 6.0/10000
    y0 = np.array([1,0])

    t,y = rk4s(f,y0,tspan,h)

    plt.xlabel("t")
    plt.ylabel("y")
    plt.plot(t,y[0],'-.',label='x(t)')
    plt.plot(t,y[1],'-',label="y(t)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
