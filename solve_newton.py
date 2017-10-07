import numpy as np
from scipy import optimize

def newton(f,df,p0,epsilon=10**-6,maxi=1000):
    """solve nonlinear eqution by newton iterative method

    :f: func
    :df: First derivative of f
    :p0: initial value
    :returns: solve

    """
    x=p0
    res=f(x)

    i=0
    while abs(res)>= epsilon:
        x1 = x - (f(x)/df(x))
        res = f(x1)
        i=i+1
        x=x1
        print("iteration %i, x=%f"%(i,x))
        if i>= maxi:
            print("maximum iteration!")
            return x

    return x


def main():
    def f(x):
        return x**3 - 3*x**2 + 4*x - 2
    def df(x):
        return 3*x**2 -6*x + 4

    x = newton(f,df,1.5)
    print("solution is ",x)

    x = optimize.root(f,1.5,jac=df)
    print("scipy root:\n",x)

if __name__ == "__main__":
    main()



