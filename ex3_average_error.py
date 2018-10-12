# -*- coding: utf-8 -*-

import numpy as np


def aver(A_):
    A = np.array(A_)
    return np.average(A)

def error(A_):
    A = np.array(A_)
    return np.std(A)/np.sqrt(len(A) - 1)


def main():
    data = [1.233, 1.242, 1.215, 1.237, 1.282, 1.243, 1.223, 1.219]
    print("the average of data is ",aver(data))
    print("the std error of data is ",error(data))


if __name__ == "__main__":
    main()
