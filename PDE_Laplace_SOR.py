# -*- coding: utf-8 -*-
"""
Created on Fri Dec  1 15:44:09 2017

@author: guilo
"""
import numpy as np
from scipy.integrate import odeint
def main():
    n=15;dpi=np.pi/n
    x = np.linspace(0,np.pi,16)
    u = np.sin(x)
    t = np.linspace(0,0.4,40)
    new_u = odeint(f,t,u)
    
    def f(t,u,a):
        x= np.zeors(len(u))
        x[0] = -2*u[0] + u[1]
        x[-1] = u[-2] - 2*u[-1]
        x[1:-1] = u[0:-4] - 2*u[1:-3] \
        +u[2:-2]
        return a*x