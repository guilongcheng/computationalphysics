import numpy as np
from scipy import linalg

def lu_decom(A,b):
    """
    LU decomposition of matrix A
    and solve Ax=b

    :A: Matrix
    :b: vector
    :returns: x, L ,U

    """
    # init
    n = len(b)
    L = np.eye(n)
    U = np.zeros((n,n))
    x = np.zeros(n)
    y = np.zeros(n)

    # decomposition A = LU

    U[0,:] = A[0,:]
    L[1:,0] = A[1:,0] / U[0,0]

    for i in range(1,n):
        for j in range(i,n):

            U[i,j] = A[i,j] - np.dot(L[i,:i],U[:i,j])

            if j != n-1:
                L[j+1,i] = (A[j+1,i] - np.dot(L[j+1,:i],U[:i,i])) / U[i,i]

    # solve Ly=b
    y[0] = b[0]

    for k in range(1,n):
        y[k] = b[k] - np.dot(L[k,:k],y[:k])

    # solve Ux=y
    x[-1] = y[-1] / U[-1,-1]

    for k in range(n-2,-1,-1):
        x[k] = (y[k] - np.dot(U[k,k+1:],x[k+1:])) / U[k,k]

    return x,L,U

def main():
    A = np.array([[2,2,3],[4,7,7],[-2,4,5]],'f8')
    b = np.array([3,1,-7],'f8')

    x,L,U = lu_decom(A,b)

    print("L is \n",L)
    print("U is \n",U)
    print("lu is \n",np.dot(L,U))
    print("solve is \n",x)

    lu,piv = linalg.lu_factor(A)
    x = linalg.lu_solve((lu,piv),b)
    print("lu is \n",lu)
    print("piv is \n",piv)
    print("solve is \n",x)




if __name__ == "__main__":
    main()

