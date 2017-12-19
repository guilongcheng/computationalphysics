import numpy as np

def gauss_jordan(A,b, eps = 1.0/(10**10)):
    n = len(b)
    dtype = A.dtype
    # print("init A is \n", A)
    for i in range(n-1):
        t = np.argmax(np.abs(A[i:,i]))+i
        if (np.abs(A[t,i]) <= 10**(-16)):
            print("matrix is singular\n" ,"index is :",t,i)
            return
        if (i != t):
            A[[i,t],:] = A[[t,i],:]
            b[i],b[t] = b[t], b[i]
            # print("t,i is : \n",t,i)
            # print("A[i,:] is \n",A[i,:])
            # print("A[t,:] is \n",A[t,:])

        m = - A[i+1:,i]/ A[i,i]
        A[i+1:,i] = 0
        A[i+1:,i+1:] = A[i+1:,i+1:] + np.outer(m,A[i,i+1:])
        b[i+1:] = b[i+1:] + m*b[i]

        # print("A[i+1:,i+1:] is \n",A[i+1:,i+1:])

    # print("final A is \n",A)
    x=np.zeros(n,dtype=dtype)

    x[n-1] = b[n-1] / A[n-1,n-1]
    for i in np.arange(n-2,-1,-1):
        x[i] = (b[i] - np.sum(x[i+1:]* A[i,i+1:])) / A[i,i]

    return x


def main():
    print("==========习题2.1.1==========")
    #Becareful the number type! if use int will not right!
    A = np.array([[1,2,1],[2,2,3],[-1,0,-3]],'f8')
    b = np.array([0,3,0],'f8')
    x = gauss_jordan(A,b)
    print("Solve is \n",x)

    x = np.linalg.solve(A,b)
    print("np solve is \n",x)
    
    print("==========习题2.1.2==========")
    #Becareful the number type! if use int will not right!
    A = np.array([[2,3,5],[3,4,8],[1,3,3]],'f8')
    b = np.array([5,6,5],'f8')
    x = gauss_jordan(A,b)
    print("Solve is \n",x)

    x = np.linalg.solve(A,b)
    print("np solve is \n",x)

if __name__ == "__main__":
    main()

