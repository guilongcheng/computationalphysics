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
    while np.sum(np.abs(res))>= epsilon:
        x1 = x - (f(x)/df(x))
        res = f(x1)
        i=i+1
        x=x1
        print("iteration %i ,"%i,"x = ",x)
        if i>= maxi:
            print("maximum iteration!")
            return x

    return x


def main():
    def f(x):
        y1 = x[1]*np.cos(x[0]*x[1]) + 1
        y2 = np.sin(x[0]*x[1]) + x[0] - x[1]
        return np.array([y1,y2],'f8')

    def df(x):
        y1 = -x[1]**2*np.sin(x[0]*x[1]) + np.cos(x[0]*x[1]) - x[1]*x[0]*np.sin(x[0]*x[1])
        y2 = x[1]*np.cos(x[0]*x[1]) + x[0]*np.cos(x[0]*x[1])
        return np.array([y1,y2],'f8')

    p0 = np.array([1,2],'f8')
    x = newton(f,df,p0)
    print("solution is ",x)

    print("scipy root \n",optimize.root(f,p0))

if __name__ == "__main__":
    main()
