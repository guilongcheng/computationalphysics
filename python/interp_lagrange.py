import numpy as np
from scipy import interpolate
import matplotlib.pyplot as plt
import time

def lagrange_1(x,y,u):
    """拉格朗日插值程序

    :x: 数据点x坐标
    :y: 数据点y坐标
    :u: 插值点x 的值
    :returns: 插值点y的值

    """
    n=len(x)
    m = len(u)
    v = np.zeros(m)

    for i in range(n):
        w = np.ones(m)
        for j in range(n):
            if j!= i:
                w = (u - x[j]) /(x[i] - x[j]) * w  # 计算拉格朗日插值基函数
        v = v + w*y[i]   #由 y(x) = \sum w_i * y_i  给出插值点函数值

    return v

def main():

    print("#====================例题3.1.5====================")
    x = np.arange(4)
    y = np.array([-5,-6,-1,16])

    u = np.arange(-0.25,3.25,0.01)
    start = time.time()
    v = lagrange_1(x,y,u)
    print("lagrange time used is %f s"%(time.time() - start))


    plt.plot(x,y,'bo',label='data')
    plt.plot(u,v,'r-',label="interpolate")
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    start = time.time()
    v = np.polyval(interpolate.lagrange(x,y),u)
    print("scipy time used is %f s"%(time.time() - start))

    plt.plot(x,y,'bo',label='data')
    plt.plot(u,v,'r-',label="np.interp")
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()

    print("===================作业3.1.8=====================")
    xi = np.linspace(0,np.pi,9)
    yi = np.cos(xi)

    x = np.linspace(0,np.pi,25)
    y = np.zeros(len(x))

    for t in range(len(x)):
        if x[t] <= xi[1]:
            y[t] =  lagrange_1(xi[0:3],yi[0:3],[x[t]])[0]
        if x[t] > xi[-2]:
            y[t] =  lagrange_1(xi[-3:],yi[-3:],[x[t]])[0]
        if x[t]>xi[1] and x[t]<= xi[-2]:
            for i in range(1,len(xi)-1):
                if xi[i-1] < x[t] <= xi[i] and abs(x[t]-xi[i-1]) <= abs(x[t]-xi[i]):
                    y[t] = lagrange_1(xi[i-2:i+1],yi[i-2:i+1],[x[t]])[0]
                if xi[i-1] < x[t] <= xi[i] and abs(x[t]-xi[i-1]) > abs(x[t]-xi[i]):
                    y[t] = lagrange_1(xi[i-1:i+2],yi[i-1:i+2],[x[t]])[0]

    plt.plot(x,y,'ro',label="lagrane 3 point")
    plt.plot(xi,yi,'bo',label='data')
    plt.plot(x,y,'r-')
    plt.grid()
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.show()






if __name__ == "__main__":
    main()



