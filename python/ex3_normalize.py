# -*- coding: utf-8 -*-
import numpy as np

def normalize(A):
    A = np.array(A)
    min_value = np.min(A)
    max_value = np.max(A)

    new_A = (A - min_value) / (max_value-min_value)
    return new_A

def main():
    A = [[3, 2, 9], [4, 7, 3], [11, 8, 5]]
    print("A is ", A)
    new_A = normalize(A)
    print("new A is ", new_A)


if __name__ == "__main__":
    main()
