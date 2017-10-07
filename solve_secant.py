import numpy as np
from scipy import optimize


def secant(f,x0,x1,eps,nmax):
    """TODO: Docstring for secant.

    :f: function
    :x0: initial x0
    :x1: initial x1
    :eps: epsilon
    :nmax: max iterations
    :returns: solution

    """
    dx = x1 -x0
    i=0

    while abs(dx) >= eps and i <= nmax:
        x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0))
        if f(x0)*f(x2) >0:
            dx = x2 -x0
            x0 = x2
        else:
            dx = x2 - x1
            x1 = x2
        i=i+1

    return x2

def main():
    def f(x):
        return x - np.exp(-x)

    x = secant(f,0,1,10**-6,100)

    print("solution is ",x)

    print("scipy root:",optimize.root(f,0.5))

if __name__ == "__main__":
    main()
