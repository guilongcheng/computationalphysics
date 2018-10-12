#!/bin/env python
import numpy as np

def TrapInt1(f,a,b,nsub):
    """TODO: Docstring for TrapInt1.

    :f: 函数
    :a: 下限
    :b: 上限
    :nsub: 区间数
    :returns: 积分结果

    """
    n=nsub+1
    h = (b-a)/nsub
    x=np.linspace(a,b,n,dtype="f8")
    y=f(x)
    T=h*(0.5*y[0]+np.sum(y[1:n-1])+0.5*y[n-1])
    return T

def SimpInt1(f,a,b,ndouble_sub):
    """复化辛普森积分公式

    :f: TODO
    :a: TODO
    :b: TODO
    :ndouble_sub: 2n个区间,n的取值
    :returns: TODO

    """
    n = 2*ndouble_sub + 1
    h = (b-a)/(n-1)
    x = np.linspace(a,b,n)
    y = f(x)
    S = (h/3)*(y[0] +4 * np.sum(y[1:n:2]) +\
               2 * np.sum(y[2:n-1:2]) + y[n-1])#注意指标！

    return S

def main():

    print("= = = = = = = = = =例题4.2.5========== ")
    def f(x):
        return 4/(1+x**2)

    print("复化梯形结果:",TrapInt1(f,0,1,16))
    print("复化辛普森结果:",SimpInt1(f,0,1,8))


    print("=====计算 integrage[1/(1+x)\sqrt(x),{x,0,1}]=====")
    def f(x):
        return 2.0/(1.0+x**2 )
    print("复化梯形结果:",TrapInt1(f,0,1,3))
    print("复化梯形结果:",TrapInt1(f,0,1,4))
    print("复化梯形结果:",TrapInt1(f,0,1,5))
    print("复化梯形结果:",TrapInt1(f,0,1,6))


if __name__ == "__main__":
    main()
