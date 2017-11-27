#!/bin/env python
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

def Eular(f,u0,a,b,h):
    """TODO: Docstring for Eular.

    :f: TODO
    :u0: TODO
    :a: TODO
    :b: TODO
    :h: TODO
    :returns: TODO

    """
    t_k = np.arange(a,b+h,h)
    n = len(t_k)
    u = np.zeros(n)
    p = np.zeros(n)
    u[0] = u0
    p[0] = u0
    for i in range(1,n):
        p[i] = p[i-1] + h*f(t_k[i-1],p[i-1])
        tmp = u[i-1] + h*f(t_k[i-1],u[i-1])
        u[i] = u[i-1] + h/2 * (f(t_k[i-1],u[i-1]) + f(t_k[i],tmp))
    return p,u


def main():
    print("==========例题5.2.1==========")
    def f(x,y):
        return -y
    p,u = Eular(f,1,0,8,0.1)
    x = np.arange(0,8+0.1,0.1)

    plt.xlabel("t")
    plt.ylabel("y(t)")
    plt.plot(x,p,'ro',label="numeric")
    plt.plot(x,np.exp(-x),'b-',label='analysis')
    plt.legend()
    plt.show()

    print("==========例题5.2.2==========")
    def f(x,y):
        return y-2*x/y
    a=0;b=1;h=0.1
    p,u = Eular(f,1,a,b,h)
    x = np.arange(a,b+h,h)

    plt.xlabel("t")
    plt.ylabel("y(t)")
    plt.plot(x,p,'ro',label="show")
    plt.plot(x,u,'go',label="hide")
    plt.plot(x,np.sqrt(1+2*x),'b-',label='analysis')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
