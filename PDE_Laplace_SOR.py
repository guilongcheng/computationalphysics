# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:44:09 2017

@author: guilongcheng
"""
import numpy as np
from scipy.integrate import odeint
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

def FTCS(kappa,u0,hx,ht,tmax):
    """TODO: Docstring for FTCS.

    :kappa: 扩散系数
    :u0: 边界条件
    :hx: 空间方向的最小划分
    :ht: 时间方向的最小划分
    :tmax: 求解最大时间
    :returns: T(x,t)

    """
    coeff = kappa*ht/(hx)**2
    if(coeff > 0.5):
        print("不满足稳定性条件！")
    nx = len(u0);nt=int(tmax/ht);
    T = np.zeros((nx,nt),'f8')
    T[:,0] = u0[:]
    T[0,:] = 0
    T[nx-1,:]=0

    for i in range(1,nt):
        T[1:nx-1,i]= T[1:nx-1,i-1] +\
            coeff * ( T[2:nx,i-1] + T[:nx-2,i-1] - 2 * T[1:nx-1,i-1])

    return T

def main():
    x = np.linspace(0,np.pi,50)
    kappa = 1
    hx = x[1]-x[0]
    ht = (hx)**2/4/kappa
    tmax = 4
    u0 = np.sin(x)

    u = FTCS(kappa,u0,hx,ht,tmax)

    x, t = np.meshgrid(np.arange(0,tmax-ht,ht),x)
    print(x.shape, t.shape,u.shape)
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.plot_surface(x,t,u)
    plt.show()

if __name__ == "__main__":
    main()
