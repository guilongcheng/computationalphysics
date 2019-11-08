# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 21:17:28 2018

@author: longcheng gui
"""
import numpy as np
import matplotlib.pyplot as plt
def spline(x, y, yp1, ypn):
# =============================================================================
#   x,y为给定的数据点，yp1和ypn为插值函数在第一个点和第n个点处的一阶导数值
#   返回数组y2为在xi 处的二阶导数值
# =============================================================================
    n = x.shape[0]
    u = np.zeros(n,'f8')
    y2 = np.zeros(n,"f8")
    if yp1 > 0.9 * 10**30:  #有yp1的值来设置边界条件，此时为自然边界条件
        y2[0] = u[0] = 0.0
    else:                   #否则有特定的一阶导数
        y2[0] = -0.5
        u[0] = (3.0/(x[1]-x[0])) * ((y[1]-y[0])/(x[1]-x[0])-yp1)

    for i in range(1,n-1): # 三对角算法的分解循环，y2和u为被分解因子的临时变量
        sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1])
        p = sig * y2[i-1] + 2.0
        y2[i] = (sig - 1.0) / p
        u[i] = (y[i+1] - y[i]) / ( x[i+1] - x[i]) - (y[i] - y[i-1])/(x[i] - x[i-1])
        u[i] = (6.0 * u[i] / (x[i+1] - x[i-1]) - sig * u[i-1]) / p

    if ypn > 0.9 * 10**30: #设定边界条件
        qn=un=0.0
    else:
        qn=0.5
        un=(3.0 / (x[n-1] - x[n-2])) * (ypn - (y[n-1] - y[n-2])/(x[n-1] - x[n-2]))

    y2[n-1] = (un-qn*u[n-2]) / (qn * y2[n-2] + 1.0)

    for k in range(n-2,-1,-1):   # 三对角算符的回代循环
        y2[k] = y2[k] * y2[k+1] + u[k]

    return y2

def splint(xdata, ydata, x, y2 = None, yp1n = None):
# =============================================================================
#     xdata,ydata 为数据点,y2为spline给出的二次导数值，x为要求的插值点坐标,yp1n为一阶导数的边界（yp1,ypn）。如果yp1n不为None，则不需要输入y2
#     返回值：在x点处的插值函数值
# =============================================================================
    n = xdata.shape[0]

    if y2 is None:
        if yp1n is not None:
            y2 = spline(xdata, ydata, yp1n[0], yp1n[1])
        else:
            print("Both of y2 and yp1n are None!!! ")
            return

    klo = 0
    khi = n-1

    y = np.zeros(x.shape[0],"f8")
    iy = 0
    for xi in x:
        # 采用二分法查找xi 位于数据点的区间。klo和khi 为xi 所在区间的上下界

        if xi > xdata[khi] or xi < xdata[klo] or iy==0:
            klo = 0
            khi = n-1
            while khi -klo > 1:
                k = (khi + klo) // 2  # // 为整除运算
                if (xdata[k] > xi):
                    khi = k
                else:
                    klo = k
            print(iy,khi,klo)

        h = xdata[khi] - xdata[klo]

        if h == 0.0:   # xdata的值不能相同
            print("Bad xdata input to routine splint")
            return 0

        # 三次样条多项式求值，即已知插值函数二次导数值来给出插值函数值
        a = (xdata[khi] - xi) / h
        b = (xi - xdata[klo]) / h
        y[iy] = a * ydata[klo] + b * ydata[khi] + ((a**3 - a) * y2[klo] +
                 (b**3 - b) * y2[khi]) * h**2 / 6.0
        iy=iy+1
    return y

def main():
    # 例题3.1.7

    xdata = np.arange(-1,1,0.1)
    ydata = 1.0 / (1.0 + xdata**2)

    x = np.arange(-0.93,0.93,0.03)
    y = splint(xdata, ydata, x, yp1n=(1.0*10**30,1.0*10**30))

    plt.plot(xdata,ydata,"bo",label="data")
    plt.plot(x,y,'r-',label="spline")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    main()
