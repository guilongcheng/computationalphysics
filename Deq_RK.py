#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt

def RungKutta(f,x_k,u0,**args):
    """

    """

    h=x_k[1]-x_k[0]
    n = len(x_k)
    y = np.zeros((len(u0),n))
    y[:,0] = u0

    for i in range(1,n):
        k1 = f(x_k[i-1], y[:,i-1],**args)
        k2 = f(x_k[i-1]+h/2, y[:,i-1]+k1*h/2,**args)
        k3 = f(x_k[i-1]+h/2, y[:,i-1]+k2*h/2,**args)
        k4 = f(x_k[i-1]+h, y[:,i-1]+h*k3,**args)
        y[:,i] = y[:,i-1] + h * (k1+2*k2+2*k3+k4)/6

    return y


def main():

    print("==========例题5.2.3==========")
    def f(x,y):
        return x*np.sqrt(y)
    a=2;b=3;u0=np.array([4]);h=0.1
    x = np.arange(a,b+h,h)
    y = RungKutta(f,x,u0)
    ye = (1+0.25*x**2)**2

    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.plot(x,y[0,:],'ro',label="numeric")
    plt.plot(x,ye,'b-',label='analysis')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
