import numpy as np
import matplotlib.pyplot as plt

def aitken(xi,fi,x):
    """逐次线性插值实现插值多项式

    :n: number of data point
    :xi: input data n point x
    :fi: input data n point y
    :x: interpolation point
    :f: interpolation value
    :df:
    :returns: None

    """
    nmax=21
    n = len(xi)
    if n> nmax:
        print("dimension too large!")
        return
    if n != len(fi):
        print("xi dimension not eq fi!")
        return

    ft = np.zeros(n,'f8')
    for i in range(n):
        ft[i] = fi[i]

    for i in range(n-1):
        # print('i:',i)
        for j in range(n-i-1):
            x1=xi[j];x2=xi[j+i+1]
            f1=ft[j];f2=ft[j+1]
            ft[j] = (x-x1)/(x2-x1)*f2 + (x-x2)/(x1-x2)*f1
            # print("j,x1,x2,f1,f2,ft:",j,x1,x2,f1,f2,ft[j])

    f=ft[0]
    df = (np.abs(f-f1)+np.abs(f-f2))/2.0
    return f,df

def main():

    xi = np.linspace(0,np.pi,9)
    fi = np.sin(xi)

    x = np.linspace(0,np.pi,21)
    f = np.zeros(len(x))
    df = np.zeros(len(x))

    for i in range(len(x)):
        f[i],df[i] = aitken(xi,fi,x[i])

    plt.plot(xi,fi,'bo',label='data')
    plt.plot(x,f,'r-',label="aitken")
    plt.plot(x,df,'y-',label="df")
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
