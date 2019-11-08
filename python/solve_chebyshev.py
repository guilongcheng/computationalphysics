from scipy import optimize
import solve_newton

def chebyshev(f,fp,fpp,p0,eps=10**-10,nmax=1000):
    """TODO: Docstring for chebyshev.

    :f: TODO
    :fp: TODO
    :fpp: TODO
    :eps: TODO
    :nmax: TODO
    :returns: TODO

    """
    x=p0
    res=f(x)
    i=0
    while abs(res)>eps:
        x1 = x - f(x)/fp(x)  - 0.5*(f(x)/fp(x))**2 * fpp(x)/fp(x)
        res = f(x1)
        i=i+1
        x = x1
        print("iteration %i, x=%f"%(i,x))
        if i>= nmax:
            print("maximum iteration!")
            return x
    return x


def chebyshev2(f,fp,p0,eps=10**-10,nmax=1000):
    """TODO: Docstring for chebyshev.

    :f: TODO
    :fp: TODO
    :eps: TODO
    :nmax: TODO
    :returns: TODO

    """
    x=p0
    res=f(x)
    i=0
    while abs(res)>eps:
        x1 = x - 0.5*f(x)/fp(x)
        x2 = x - f(x)/fp(x1)
        res = f(x2)
        i=i+1
        x = x2
        print("iteration %i, x=%f"%(i,x))
        if i>= nmax:
            print("maximum iteration!")
            return x
    return x

def main():
    def f(x):
        return x**3 - 2
    def fp(x):
        return 3*x**2
    def fpp(x):
        return 6*x

    print("chebyshev solution is ",chebyshev(f,fp,fpp,1.0))
    print("chebyshev2 solution is ",chebyshev2(f,fp,1.0))
    print("newton solution is ",solve_newton.newton(f,fp,1.0,epsilon=10**-10))

    print("scipy root:\n",optimize.root(f,1,jac=fp))

if __name__ == "__main__":
    main()
