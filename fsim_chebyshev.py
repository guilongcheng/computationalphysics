import numpy as np
import matplotlib.pyplot as plt
try :
    from mtime import mtime
except ImportError:
    print("Package mtime not find!")


def chebft(a, b, n, func):
    """
    a,b 是 x 的取值范围[a,b]，
    n 为 系数的个数
    func 为给定函数
    """
    k = np.arange(n)
    y = np.cos(np.pi*(k+0.5)/n) # cos(pi *( k - 1/2) /N )
    f = func( y * (b-a)/2 + (b+a)/2)  # 将 y 转化到[a,b]上

    j = np.arange(n)

    c = np.dot(np.cos(np.pi * j.reshape(n,1) * (k.reshape(1,n)+0.5) /n ), f)
    c = c*2.0/n
    return c

def chebev(a, b, c, m, x):
    """
    a,b 是x的取值范围
    c为切比雪夫系数数组
    m为近似的项数
    x为需求取的坐标值
    """
    if np.any((x-a)*(x-b) > 0):
        print("x 不在区间(a,b)中")
        return 0

    y = (2.0 * x -a -b) / (b -a)   # 变量替换到(-1,1)区间
    y2 = 2.0 * y

    d=np.zeros(len(x),'f8')
    dd=np.zeros(len(x),'f8')
    for j in range(m-1,0,-1):
        sv = d
        d = y2 * d - dd + c[j]
        dd = sv

    return y*d - dd + 0.5*c[0]

def main():
    def f1(x):
        return np.sin(x**0.5)/ x**0.5

    a=0;b=2*np.pi;n=m=3
    x = np.linspace(a+0.01,b-0.01,30)

    c = chebft(a,b,n,f1)

#    t = mtime()

#    t.start("sim")
    sim_data = chebev(a,b,c,m,x)
#    t.end("sim")

#    t.start("real")
    real_data = f1(x)
#    t.end("real")

    print("abs(Delta f)  is ",np.sum(np.abs(real_data-sim_data) ))

    plt.plot(x,real_data,'-')
    plt.plot(x,sim_data,'o')
    plt.show()


if __name__ == "__main__":
    main()




