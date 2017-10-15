import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import time

def lagrange_1(x,y,u):
    """拉格朗日插值程序

    :x: TODO
    :y: TODO
    :u: TODO
    :returns: TODO

    """
    n=len(x)
    m = len(u)
    v = np.zeros(m)

    for i in range(n):
        w = np.ones(m)
        for j in range(n):
            if j!= i:
                w = (u - x[j]) /(x[i] - x[j]) * w
        v = v + w*y[i]

    return v

def main():
    x = np.arange(4)
    y = np.array([-5,-6,-1,16])

    u = np.arange(-0.25,3.25,0.01)
    start = time.time()
    v = lagrange_1(x,y,u)
    print("lagrange time used is %f s"%(time.time() - start))


    plt.plot(x,y,'bo',label='data')
    plt.plot(u,v,'r-',label="interpolate")
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    start = time.time()
    v = np.polyval(interpolate.lagrange(x,y),u)
    print("scipy time used is %f s"%(time.time() - start))

    plt.plot(x,y,'bo',label='data')
    plt.plot(u,v,'r-',label="np.interp")
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()



