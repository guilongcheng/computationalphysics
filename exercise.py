import numpy as np
from scipy import integrate
import time
from integ_Romberg import Roberg

def main():
    pi = np.pi
    theta0 = np.array([pi/2**i for i in range(1,8)])
    k,flag=eval(input("choice theta0 and method('sp or glc')\n"))
    print("theta0 is ",theta0[k])

    l=0.5
    g=9.8

    def func(x):
        return 1/np.sqrt(np.cos(x) - np.cos(theta0[k]))

    res = 1
    i = 1
    r0 = theta0[k]*(1-1/10)
    maxi=1000
    integ0=0
    time0 = time.time()
    while res > 10**(-3) and i<maxi:
        if "sp" in flag:
            # integ=integrate.romberg(func,0,r0,tol=10**(-6),divmax=20)
            integ,tmp=integrate.quad(func,0,r0)
        if "glc" in flag:
            integ=Roberg(func,0,r0,10**(-6))
        i=i+1
        r0 = theta0[k]*(1-1/10**i)
        res = abs(integ - integ0)
        integ0 = integ

    print("T is ",4*np.sqrt(l/(2*g))*integ0)
    print("small-angle approximation is ", 2*pi * np.sqrt(l/g))
    print("all time used is ",time.time()-time0)

if __name__ == "__main__":
    main()


