import numpy as np

def gsdl(A,b,p0=None,delta=10**-5,max1=100):
    """Gauss-Seidel iterative method

    :A: Matrix
    :b: vector
    :delta: accuracy
    :max: Maximum iterations
    :returns: solve

    Add the Pivoting code by youself.
    """
    n = len(b)

    if p0 is None:
        x = np.zeros(n)
    else:
        x[:] = p0[:]

    y = np.zeros(n)

    for iters in range(max1):
        y[:] = x[:]

        for i in range(n):
            z = b[i]
            for j in range(n):
                if j != i:
                    z = z- A[i,j] * x[j]

            if np.abs(A[i,i]) < 10**(-10) or iters == max1 -1 :
                print("Fail!")
                return y

            z = z / A[i,i]
            x[i] = z

        if np.sum(np.abs(y-x)) < delta:
            return y

def main():

    A = np.array([[5,1,0],[1,5,1],[0,1,5]],'f8')
    b = np.array([4,-3,4],'f8')

    x = gsdl(A,b)

    print("solve is \n",x)

    x = np.linalg.solve(A,b)
    print("np solve is \n",x)

if __name__ == "__main__":
    main()
