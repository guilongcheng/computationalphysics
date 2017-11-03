import numpy as np
import sympy as sp
from numpy import sin,cos,exp


def filon(f,a,b,k,eps):
    """filon 计算震荡积分

    :f: TODO
    :a: TODO
    :b: TODO
    :eps: TODO
    :returns: TODO

    """

    nsub = 2
    C=0
    S=0
    tol1 = 1
    tol2 = 1
    while tol1>eps and tol2 >eps:
        C0 = C
        S0 = S

        nsub = 2*nsub
        n = nsub+1
        h = (b-a) / nsub
        x = np.linspace(a,b,n)
        theta = k*h

        alpha = 1/theta + sin(2*theta)/(2*theta**2) - (2*sin(theta)**2)/theta**3
        beta = 2*( (1+cos(theta)**2)/theta**2 - sin(2*theta)/theta**3 )
        gamma = 4*(sin(theta)/theta**3 - cos(theta)/theta**2)

        yi = f(x)

        C_even = np.sum(yi[0:-1:2]*cos(k*x[0:-1:2])) - 0.5*(yi[0]*cos(k*x[0]) + yi[-1]*cos(k*x[-1]))
        C_odd = np.sum(yi[1:-2:2]*cos(k*x[1:-2:2]))

        S_even = np.sum(yi[0:-1:2]*sin(k*x[0:-1:2])) - 0.5*(yi[0]*sin(k*x[0]) + yi[-1]*sin(k*x[-1]))
        S_odd = np.sum(yi[1:-2:2]*sin(k*x[1:-2:2]))

        C = h*(alpha*(yi[-1]*sin(k*x[-1]) - yi[0]*sin(k*x[0])) + beta*C_even + gamma*C_odd )
        S = h*(alpha*(yi[0]*cos(k*x[0]) - yi[-1]*cos(k*x[-1])) + beta*S_even + gamma*S_odd )
        tol1 = abs(C-C0)
        tol2 = abs(S-S0)

    print("n is ",nsub)
    return S,C


def main():
    print("===========例题4.2.9==========")
    def f(x):
        return exp(-x)

    print("answer is ",filon(f,0,10,100,10**-6))


if __name__ == "__main__":
    main()
