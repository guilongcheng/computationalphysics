import numpy as np
from scipy import optimize
import matplotlib.pyplot as plt
# from matplotlib.font_manager import FontProperties

def linearfit(x,y):
    """fit y = ax+b

    :x: data coordinate
    :y: data
    :returns: a,b

    """
    n     = len(x)
    xsum  = np.sum(x,0)
    x2sum = np.sum(x**2,0)
    ysum  = np.sum(y,0)
    xysum = np.dot(x,y)

    b = (ysum*x2sum-xsum*xysum)/(n*x2sum-xsum**2)
    a = (n*xysum-xsum*ysum)/(n*x2sum - xsum**2)

    return a,b

def f(x,a,b):
    """定义线性函数"""
    return a*x + b

def main():
    print("==========例题3.2.1==========")
    x = np.array([0.5, 1.2, 2.1, 2.9, 3.6, 4.5, 5.7])
    y = np.array([2.81, 3.24, 3.80, 4.30, 4.73, 5.29, 6.03])

    a,b = linearfit(x,y)
    print("linear parameter is ",a,b)

    plt.plot(x,y,'o',label='data')
    plt.plot(x,f(x,a,b),'-',label="fit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    popt,pcov = optimize.curve_fit(f,x,y)
    print("curve_fit linear parameter is ",popt[0],popt[1])

    plt.plot(x,y,'o',label='data')
    plt.plot(x,f(x,*popt),'-',label="curve_fit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    print("==========作业3.2.4==========")
    x = np.array([0.3,0.5,0.6,0.7,0.9])
    sinx = np.sin(x)
    y = np.array([1.37731,1.48766,1.53879,1.58653,1.67])

    a,b = linearfit(sinx,y)
    print("linear parameter is ",a,b)

    plt.plot(x,y,'o',label='data')
    plt.plot(x,f(sinx,a,b),'-',label="fit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()


