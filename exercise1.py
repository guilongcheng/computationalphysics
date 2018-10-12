import numpy as np
from scipy import integrate
import time
from integ_Romberg import Roberg

def func(x,theta0):
    a = np.cos(theta0)
    if np.min(np.abs(np.sqrt(np.cos(x)-a))) < 10**(-16):
        result = 1/(1-a) + x**2/(2*(a-1)**2) + (-a-5)*x**4/(24*(a-1)**3)
    else:
        result=1/np.sqrt(np.cos(x) - a )
    return result

def main():
    pi = np.pi
    #theta0 = np.array([pi/2**i for i in range(1,8)])
    theta0=eval(input("input theta0 \n"))
    print("theta0 is  = %.15f  %f*pi"%(theta0,theta0/np.pi))

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

    while res > eps and i<maxi:
        integ=Roberg(func,r0,r1,10**(-7),theta0=theta0)
        i=i+1
        r0 = r1
        r1 = theta0*(1-1/10**i)
        res = abs(integ)
        suminteg = suminteg + integ

    print("solved by RK:")
    print("T is ",4*np.sqrt(l/(2*g))*suminteg)
    print("all time used is ",time.time()-time0)

    time0 = time.time()
    suminteg,tmp=integrate.quad(func,0,theta0,args=(theta0,))
    print("solved by sp:")
    print("T is ",4*np.sqrt(l/(2*g))*suminteg)
    print("all time used is ",time.time()-time0)

    print("small-angle approximation is ", 2*pi * np.sqrt(l/g))
if __name__ == "__main__":
    main()


