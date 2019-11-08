import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


def lsqfit(f,dyda,p0,x,y,sigma=None,nmax=1000,delta=0.1):
    """TODO: Docstring for lsqfit.

    :f: fit function
    :x: x data
    :y: y data
    :sigma: error of y
    :nmax: max iter
    :delta: stop parameters
    :returns: fit parameters

    """

    def chi2(f,p,x,y,sigma):
        return np.sum((f(x,*p) - y)**2 / sigma**2,0)

    lamda = 0.001

    ma = len(p0)
    ndata = len(x)

    if sigma is None:
        sigma = np.ones(ndata)

    res = chi2(f,p0,x,y,sigma)
    beta = np.zeros(ma)
    alpha = np.zeros((ma,ma))

    a=p0

    iters = 0

    while res > delta or res <0 :
        beta = np.sum(((y-f(x,*a))/sigma**2).reshape((ndata,1)) * dyda(x,*a),0)
        tmp = dyda(x,*a)/sigma.reshape((ndata,1))
        alpha = np.einsum("ij,ik",tmp,tmp)
        for i in range(ma):
            alpha[i][i] = alpha[i][i]*(1 + lamda)

        da = np.linalg.solve(alpha,beta)

        tmp1 = chi2(f,a,x,y,sigma)
        tmp2 = chi2(f,a+da,x,y,sigma)
        print("tmp2,tmp1",tmp2,tmp1)
        if tmp2 >= tmp1:
            lamda = 10*lamda
            res = tmp1 - tmp2
        if tmp2 < tmp1 :
            lamda = 0.1*lamda
            res = tmp1 - tmp2
            a = a + da
        iters=iters+1
        if (i >= nmax):
            print("not conver")
            return 0

    print("all iter is ",iters)


    return a


def main():

    x = np.arange(0,10*np.pi/3,0.1)
    ndate = len(x)
    nconf=100
    y = np.zeros((ndate,nconf))
    for i in range(ndate):
        y[i,:] = np.random.normal(loc=np.sin(0.3*x[i])+4, scale=0.1,size=nconf)
    sigma = np.std(y,1) / (nconf-1)**0.5
    y_ave = np.average(y,1)

    def func(x,a,b):
        return np.sin(a*x) + b

    def dyda(x,a,b):
        return np.array([x*np.cos(a*x),np.ones(len(x))]).T

    params = lsqfit(func,dyda,(0.2,0.2),x,y_ave,sigma=sigma,delta=0.01)

    plt.errorbar(x,y_ave,sigma,fmt='o',label='data')
    plt.plot(x,func(x,*params),'-',label="fit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    print("lsqfit params:",params)

    popt,pcov = curve_fit(func,x,y_ave,p0=(0.2,0.2),sigma=sigma)

    plt.errorbar(x,y_ave,sigma,fmt='o',label='data')
    plt.plot(x,func(x,*popt),'-',label="fit")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    print("curve_fit:",popt)

if __name__ == "__main__":
    main()






