import numpy as np
from scipy import integrate
import time
from integ_Romberg import Roberg

def func(x,theta0):
    try:
        result=1/np.sqrt(np.cos(x)-np.cos(theta0))
    except RuntimeWarning:
        a = np.cos(theta0)
        result = 1/(1-a) + x**2/(2*(a-1)**2) + (-a-5)*x**4/(24*(a-1)**3)
    return result

def main():
    pi = np.pi
    #theta0 = np.array([pi/2**i for i in range(1,8)])
    theta0,flag=eval(input("input theta0 and method('sp or glc')\n"))
    print("theta0 is  = %.15f  %f*pi"%(theta0,theta0/np.pi))

    if "glc" in flag:
        eps=10**(-6)
    if "sp" in flag:
        eps=10**(-6)

    l=0.5
    g=9.8

    res = 1
    i = 1
    r0 = 0
    r1 = theta0*(1-1/10)
    maxi=1000
    suminteg=0
    time0 = time.time()

    if "sp" in flag:
        # integ=integrate.romberg(func,0,r0,tol=10**(-6),divmax=20)
        suminteg,tmp=integrate.quad(func,0,theta0,args=(theta0,))
    if "glc" in flag:
        while res > eps and i<maxi:
            integ=Roberg(func,r0,r1,10**(-7),theta0=theta0)
            i=i+1
            r0 = r1
            r1 = theta0*(1-1/10**i)
            res = abs(integ)
            suminteg = suminteg + integ

    print("T is ",4*np.sqrt(l/(2*g))*suminteg)
    print("small-angle approximation is ", 2*pi * np.sqrt(l/g))
    print("all time used is ",time.time()-time0)

if __name__ == "__main__":
    main()


