"""该程序用来求解谐振子系统取k=1
采用最简单的预修正方法"""
import numpy as np
import matplotlib.pyplot as plt
from math import ceil

from Deq_Eular import Eular


def Pre_correct(f, u0, a, b, h=None, n=None,flag="Eular"):
    if n is None:
        n = ceil((b-a) / h)
        x = np.arange(a, b, h)
    if h is None:
        h = (b -a) / n
        x = np.linspace(a, b, n)
    m = len(u0)
    y = np.zeros((n,m))
    y[0] = u0 
    for i in range(1,n):
        y[i] = y[i-1] + f(x[i-1],y[i-1])* h 
        if flag is not "Eular": 
            y[i] = y[i-1] + h/2 * (f(x[i-1],y[i-1]) + f(x[i],y[i])) 
    
    return x,y

def main():
    "计算简谐振动k=1 m=1"
    def f(x,y):
        return np.array([y[1], -y[0]]) 
    a=0;b=4*np.pi
    u0=np.array([0,1])
    x,y1 = Pre_correct(f,u0,a,b,h=0.02*np.pi)
    x,y2 = Pre_correct(f,u0,a,b,h=0.02*np.pi,flag="Pre")

    plt.plot(x,y1[:,0],"bo",label="y(x)-Eular")
    plt.plot(x,y1[:,1],"bo",label="y'(x)-Eular")
    plt.plot(x,y2[:,0],"rs",label="y(x)-pre")
    plt.plot(x,y2[:,1],"rs",label="y'(x)-pre")
    plt.plot(x,np.sin(x),"k-",label="y(x)-ana")
    plt.plot(x,np.cos(x),"k-",label="y'(x)-ana")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()    