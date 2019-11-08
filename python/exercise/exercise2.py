# -*- coding: utf-8 -*-
"""
Created on Fri Nov 24 16:51:13 2017

@author: guilongcheng
"""

import numpy as np
from scipy.integrate import simps
import matplotlib.pyplot as plt
from Deq_RK import RungKutta

def potential(x,y,k):
    return np.array([y[1],-k**2*y[0]])

def main():
    #测试
    print("==========求解一维定态薛定谔方程==========")
    def func(x,y,V,k):
        return np.array([y[1],-(k-V(x))*y[0]])
    x = np.linspace(0,1,1000)
    k = 1
    m=3
    eigenV = np.zeros(m)
    tol = 10**(-8)
    for i in range(m):
            dk = 1/20
            k = k + dk
            u0 = np.array([0,0.001])
            phi = RungKutta(func,x,u0,k=k)
            oldphi = phi[0,-1]
            dphi = oldphi
            while abs(dphi) > tol:
                k = k+dk
                phi = RungKutta(func,x,u0,k=k)
                dphi = phi[0,-1]
                if dphi*oldphi < 0:
                    k = k - dk
                    dk = dk/2
            eigenV[i] = k
            norm = simps(phi[0,:]**2)**0.5
            print(norm)
            plt.plot(x,phi[0,:]/norm,label="value[%i]=%f"%(i,eigenV[i]))
            plt.legend()
            plt.show()

if __name__ == "__main__":
    main()





