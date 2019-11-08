import numpy as np
from scipy import optimize

def iter(f,x0,eps=10**-6,nmax = 100):

    i = 0
    x = f(x0)

    while abs(x-x0) >= eps:
        x0=x;x=f(x0);i=i+1
        print("iteration %i : x = %f"%(i,x))

    return x

def main():
    def f(x):
        return 1/2 + np.sin(x)

    x = iter(f,1)
    print("solution is ",x)

    def f2(x):
        return 1/2 + np.sin(x) - x
    print("scipy root\n",optimize.root(f2,1))

if __name__ == "__main__":
    main()


