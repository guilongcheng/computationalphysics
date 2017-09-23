import numpy as np

def gauss_jordan(A,b, eps = 1.0/(10**10)):
    n = len(b)
    dtype = A.dtype
    for i in range(n-1):
        t = np.argmax(np.abs(A[i:,i]))
        if (np.abs(A[t,i]) <= 10**(-16)):
            print("matrix is singular\n" ,"index is :",t,i)
            return
        if (i != t):
            A[[i,t],:] = A[[t,i],:]
            b[i],b[t] = b[t], b[i]
            print("t,i is :",t,i)
            print("A[i,:] is ",A[i,:])
            print("A[t,:] is ",A[t,:])

        m = - A[i+1:,i]/ A[i,i]
        A[i+1:,i] = 0
        A[i+1:,i+1:] = A[i+1:,i+1:] + np.outer(m,A[i,i+1:])
        b[i+1:] = b[i+1:] + m*b[i]

        # for j in range(i+1,n):
            # m = - A[j,i]/ A[i,i]
            # A[j,i] = 0
            # A[j,i+1:] = A[j,i+1:] + m * A[i,i+1:]
            # b[j] = b[j] + m*b[i]

        print("A[i+1:,i+1:] is ",A[i+1:,i+1:])

    x=np.zeros(n,dtype=dtype)

    x[n-1] = b[n-1] / A[n-1,n-1]
    for i in np.arange(n-2,-1,-1):
        x[i] = (b[i] - np.sum(x[i+1:]* A[i,i+1:])) / A[i,i]

    return x

