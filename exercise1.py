import numpy as np
from scipy import integrate
import time
from integ_Romberg import Roberg

def main():
    pi = np.pi
    theta0 = np.array([pi/2**i for i in range(1,8)])
    k,flag=eval(input("choice theta0 and method('sp or glc')\n"))
    print("theta0 is pi/%i = %.15f"%(2**(k+1),theta0[k]))

    if "glc" in flag:
        eps=10**(-3)
    if "sp" in flag:
        eps=10**(-6)

    l=0.5
    g=9.8

    def func(x,theta0):
        return 1/np.sqrt(np.cos(x) - np.cos(theta0))

    res = 1
    i = 1
    r0 = theta0[k]*(1-1/10)
    maxi=1000
    integ0=0
    time0 = time.time()
    while res > eps and i<maxi:
        if "sp" in flag:
            # integ=integrate.romberg(func,0,r0,tol=10**(-6),divmax=20)
            integ,tmp=integrate.quad(func,0,r0)
        if "glc" in flag:
            integ=Roberg(func,0,r0,10**(-6),theta0=theta0[k])
        i=i+1
        r0 = theta0[k]*(1-1/10**i)
        print("theta0' is ",r0)
        res = abs(integ - integ0)
        integ0 = integ

    print("T is ",4*np.sqrt(l/(2*g))*integ0)
    print("small-angle approximation is ", 2*pi * np.sqrt(l/g))
    print("all time used is ",time.time()-time0)

if __name__ == "__main__":
    main()


