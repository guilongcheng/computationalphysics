# -*- coding: utf-8 -*-

import numpy as np


def dot(A_,B_):
    A = np.array(A_)
    B = np.array(B_)
    if A.shape[1] != B.shape[0]:
        print("Error: the number of columns of A not eq the number of rows of B ")
        return 0
    C = np.zeros((A.shape[0],B.shape[1]))
    for i in range(A.shape[0]):
        for j in range(B.shape[1]):
            for k in range(A.shape[1]):
                C[i,j] = C[i,j] + A[i,k] * B[k,j]
    return C



def main():
    A = [[3.0, 4, 5], [2, 7, 9], [1, 6, 8]]
    B = [[2.0, 9, 1], [1, 5, 2], [3, 7, 4]]
    C = dot(A,B)
    print("loop A dot B as :",C)

    C = np.dot(np.array(A),np.array(B))
    print("the result of np.dot as :",C)

if __name__ == "__main__":
    main()
