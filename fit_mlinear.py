import numpy as np
# from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from linalg_gauss_jordan import gauss_jordan

def multifit(x,y,m):
    """TODO: Docstring for multifit.

    :arg1: TODO
    :returns: TODO

    """
    n = len(x)
    if m > n :
        print("error: m >n")

    S = np.zeros((m+1,m+1))
    T = np.zeros(m+1)
    tmp = np.zeros(2*m+1)

    for i in range(2*m+1):
        tmp[i] = np.sum(x**i,0)

    for i in range(m+1):
        S[i,:] = tmp[i:m+1+i]

    for i in range(m+1):
        T[i] = np.sum(x**i*y,0)

    a = gauss_jordan(S,T)

    return a


def main():

    x = np.linspace(-1,1,11)
    y = np.exp(x)


    m=5
    a = multifit(x,y,m)
    a = a[::-1]

    plt.plot(x,y,'o',label='data')
    plt.plot(x,np.polyval(a,x),'-',label="fit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    x = np.linspace(-5,2,50)
    y = np.exp(x)

    m=5
    a = multifit(x,y,m)
    a = a[::-1]

    plt.plot(x,y,'o',label='data')
    plt.plot(x,np.polyval(a,x),'-',label="fit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()
