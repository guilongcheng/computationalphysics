#!/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from decimal import *
import math
import time

def accurace(pi,flag=""):
    print("pi value is ", pi )
    if "Decimal" in flag:
        print("accurance is ", (Decimal(pi) - Decimal(math.pi))/Decimal(math.pi) *Decimal(100))
    else:
        print("accurance is ", (pi - math.pi)/math.pi *100)

def Montecarlo(N):
    data = np.array([np.random.random(N), np.random.random(N)])
    Rdata = data[0,:]**2 + data[1,:]**2
    N_r = len(np.where(Rdata <= 1.0)[0])
    pi = N_r / N * 4
    accurace(pi)

def Series1(N):
    a_s1 = 1 / (np.arange(1,N)**2)
    pi =(np.sum(a_s1)*6)**0.5
    accurace(pi)

def BaileyBorweinPlouffe(N):
    list1 = 8*np.arange(N)
    list16 = np.array([16.0**i for i in np.arange(N)])
    pi = np.sum(1.0/list16*(4/(list1+1) - 2/(list1+4) - 1/(list1+5) - 1/(list1+6)))
    accurace(pi)


def BaileyBorweinPlouffeH(N):
    N = int(N)
    list1 = 8*np.arange(N)
    list16 = np.array([16**i for i in range(N)])
    list2 = [Decimal(1/16**i)*((Decimal(4)/(8*i+1)) - (Decimal(2)/(8*i+4)) - (Decimal(1)/(8*i+5)) - (Decimal(1)/(8*i+6))) for i in range(N)]
    pi = np.sum((Decimal(1)/list16)*((Decimal(4)/(list1+1)) - (Decimal(2)/(list1+4)) - (Decimal(1)/(list1+5)) - (Decimal(1)/(list1+6))))
    accurace(pi,"Decimal")
    pi = np.sum(list2)
    accurace(pi,"Decimal")

def main():
    N = eval(input("please input N\n"))
    t_i = time.time()
    Montecarlo(N)
    t_f = time.time()
    print("all time used is ", t_f -t_i)

    t_i = time.time()
    Series1(N)
    t_f = time.time()
    print("all time used is ", t_f -t_i)

    t_i = time.time()
    BaileyBorweinPlouffe(N/100)
    t_f = time.time()
    print("all time used is ", t_f -t_i)

    t_i = time.time()
    BaileyBorweinPlouffeH(N/100)
    t_f = time.time()
    print("all time used is ", t_f -t_i)
if __name__ == "__main__":
    main()
