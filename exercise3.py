# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:51:13 2017

@author: guilongcheng
"""

import numpy as np
import matplotlib.pyplot as plt
from Deq_RK import RungKutta
from scipy.integrate import odeint

def main():
    #测试
    print("==========求解高维线性方程组==========")
    def func(t,y):
        return np.array([y[1],y[0]-y[2]-(3*y[1])**2+y[3]**3+6*y[4]+2*t,\
                         y[3],y[4],y[4]-y[1]+np.exp(y[0])-t])
    u0=np.array([2,-4,-2,7,6])
    x = np.linspace(1,1.6,1000)
    y = RungKutta(func,x,u0)
	
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.plot(x,y[0,:],'r-.',label="x(t)")
    plt.plot(x,y[2,:],'b-.',label="y(t)")
    plt.legend()
    plt.show()
   
    def func2(y,t):
        return np.array([y[1],y[0]-y[2]-(3*y[1])**2+y[3]**3+6*y[4]+2*t,\
                         y[3],y[4],y[4]-y[1]+np.exp(y[0])-t])
    sol = odeint(func2,u0,x)
    plt.xlabel("x")
    plt.ylabel("y(x)")
    plt.plot(x,sol[:,0],'r-',label="x(t)")
    plt.plot(x,sol[:,2],'b-',label="y(t)")
    plt.legend()
    plt.show()
    
if __name__ == "__main__":
    main()





