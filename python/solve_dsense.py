import numpy as np
from scipy import optimize

def dsense(myfun,x0,eps=10**-6):
    """TODO: Docstring for dsense.

    :myfun:
    :x0: TODO
    :eps: TODO
    :returns: TODO

    """
    x=x0
    f,f2,g,g2 = myfun(x)
    lamda = f2/ g2
    x = x - lamda*g
    i=0

    while np.sum(np.abs(x-x0))>=eps:
        x0 = x ; i = i+1
        f,f2,g,g2 = myfun(x)
        lamda = f2/g2
        x = x - lamda*g
        if i > 100000:
            print("don't convergence!")
            return x
    return x

def main():
    def myfun(x):
        fun1 = x[0]**2 + 4* x[1]**2 - 9
        fun2 = 18*x[1] - 14*x[0]**2 + 45
        f = np.array([fun1,fun2],'f8')
        f2 = fun1**2+fun2**2

        gf1 = 4*x[0]*fun1 -56*x[0]*fun2
        gf2 = 16*x[1]*fun1 +36*fun2
        g = np.array([gf1,gf2],'f8')
        g2 = gf1**2+gf2**2
        return f,f2,g,g2

    x0=np.array([1,1],'f8')

    x= dsense(myfun,x0,eps=10**-13)
    print("solution is ",x)

    def fmin(x):
        fun1 = x[0]**2 + 4* x[1]**2 - 9
        fun2 = 18*x[1] - 14*x[0]**2 + 45
        return fun1**2 + fun2**2

    def fjac(x):
        fun1 = x[0]**2 + 4* x[1]**2 - 9
        fun2 = 18*x[1] - 14*x[0]**2 + 45
        gf1 = 4*x[0]*fun1 -56*x[0]*fun2
        gf2 = 16*x[1]*fun1 +36*fun2
        g = np.array([gf1,gf2],'f8')
        return g

    print("scipy root:\n",optimize.minimize(fmin,x0,jac=fjac))

if __name__ == "__main__":
    main()

