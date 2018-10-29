import numpy as np
import mtime
from scipy import linalg

def tridag(a,b,c,f):
    """TODO: Docstring for tridag.

    :arg1: TODO
    :returns: TODO

    """
    n = len(f)
    x = np.zeros(n,'f8')
    y = np.zeros(n,'f8')
    d = np.zeros(n,'f8')
    u = np.zeros(n,'f8')

    d[0] = b[0]

    for i in range(n-1):
        u[i] = c[i] / d[i]
        d[i+1] = b[i+1] - a[i+1] * u[i]

    y[0] = f[0] / d[0]

    for i in range(1,n):
        y[i] = (f[i] - a[i]*y[i-1])/d[i]

    x[n-1] = y[n-1]

    for i in np.arange(n-2,-1,-1):
        x[i] = y[i] - u[i] * x[i+1]

    return x


def main():
    a = -np.ones(100,'f8')
    b = 2*np.ones(100,'f8')
    c = -np.ones(100,'f8')
    f = np.arange(1,101,dtype='f8')

    t = mtime.mtime()
    t.start("tridag")
    x = tridag(a,b,c,f)
    t.end("tridag")
    print("sovle is \n",x)

# one way to construct matrix
    d = np.zeros((100,100))
    for i in range(100):
        d[i,i] = b[i]
        if i <99:
            d[i+1,i] = a[i+1]
            d[i,i+1] = c[i]

    t.start("numpy solve")
    x = np.linalg.solve(d,f)
    t.end("numpy solve")

    t.start("scipy solve")
    x = linalg.solve(d,f)
    t.end("scipy solve")
    print("solve is \n",x)

# other way to construct matrix
    d = np.diagflat(b) + np.diagflat(a[:-1],1) + np.diagflat(c[:-1],-1)
    ab = np.array([a,b,c])

    t.start("scipy solve_banded")
    x = linalg.solve_banded((1,1),ab,f)
    t.end("scipy solve_banded")

    print("solve is \n",x)

if __name__ == "__main__":
    main()
