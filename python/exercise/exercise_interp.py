import numpy as np
import matplotlib.pyplot as plt
from mtime import mtime



def lagrange3(x,y,u):
    v = np.zeros(len(u),u.dtype)

    def calcv(u,x,y):
        deltax1 = x[1] - x[0]
        deltax2 = x[2] - x[1]
        deltax3 = x[2] - x[0]
        deltau1 = u - x[0]
        deltau2 = u - x[1]
        deltau3 = u - x[2]
        return deltau2*deltau3 /(deltax3*deltax1) * y[0] \
                - deltau1*deltau3 / (deltax1*deltax2) * y[1] \
                + deltau1*deltau2 / (deltax3*deltax2) * y[2]


    for i in range(len(u)):
        if u[i] < x[1]:
            v[i] = calcv(u[i],x[:3],y[:3])
        elif u[i] > x[-2]:
            v[i] = calcv(u[i],x[-3:],y[-3:])
        else:
            for j in range(2,len(x)-1):
                if u[i] < x[j]:
                    if abs(u[i]-x[j-1]) <= abs(u[i]-x[j]):
                        v[i] = calcv(u[i],x[j-2:j+1],y[j-2:j+1])
                    else:
                        v[i] = calcv(u[i],x[j-1:j+2],y[j-1:j+2])
                    break
    return v
def main():

    print("#====================例题3.1.3====================")
    x = np.linspace(0,2*np.pi,9) # 数据点x坐标
    y = np.sin(x)                # 数据点y坐标


    u = np.linspace(0,2*np.pi,21) # 返回新插值函数x坐标
    t = mtime()                    # 记录时间
    t.start("lagrange")
    v = lagrange3(x,y,u)           # 插值函数计算
    t.end("lagrange")

    #  画图
    plt.plot(x,y,'bo',label='data')        # 数据点
    plt.plot(u,v,'r-',label="interpolate") # 插值函数
    plt.grid()                             # 画网格
    plt.xlabel("x")                        # x坐标标签
    plt.ylabel("y")                        # y坐标标签
    plt.legend()                           # 图例
    plt.show()                             # 显示图片

if __name__ == "__main__":
    main()
