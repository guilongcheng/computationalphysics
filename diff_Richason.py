#!/bin/env python
import numpy as np

def Richason(f,x0,n,h):
    """理查森外推法

    :f: TODO
    :x0: TODO
    :n: TODO
    :h: TODO
    :returns: TODO

    """
    G = np.zeros(n)
    for i in range(n):
        y1 = f(x0+h/2**(i+1))
        y2 = f(x0-h/2**(i+1))
        G[i] = 2**(i)*(y1-y2)/h

    G1=np.copy(G)

    for i in range(n-1):
        for j in range(i+1,n):
            G1[j] = (G[j] - 0.5**(2*(i+1))*G[j-1]) / (1 - 0.5**(2*(i+1)))

        G[:] = G1[:]

    df = G[-1]
    return df

def main():

    def f(x):
        return 2**x
    print("y=2**x 在x=1处的导数为:%.15f, 精确值:%.15f"%(Richason(f,1,8,0.1),2*np.log(2)))

if __name__ == "__main__":
    main()
