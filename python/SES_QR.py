# -*- coding: utf-8 -*-

import numpy as np

from numpy.linalg import qr as npqr



def OQR(A,nmax=40):
    N = A.shape[0]
    QM=np.eye(N)

    T = A
    for i in range(nmax):
        q,r = npqr(T-np.eye(N)*T[-1,-1])
        T = np.dot(r,q)+np.eye(N)*T[-1,-1]

    return T






def main():

    np.set_printoptions(precision=4,suppress=True)
    A = np.array([[3,17,-37,18,-40],[1,0,0,0,0],[0,1,0,0,0],[0,0,1,0,0],[0,0,0,1,0]])
    A = np.array([[17,24,-1,8,15],[23,5,7,14,16],[4,6,13,20,22],[10,12,19,21,3],[11,18,25,2,9]])


    q,r= npqr(A)

    print("q is ",q)
    print("r is ",r)
    h,tau= npqr(A,mode="raw")
    print("h is ",h)
    print("tau is ",tau)
    print("h.A is ",np.dot(np.transpose(h),A))

    T = OQR(A,100)
    print("T is ",T)

if __name__ == "__main__":
    main()

