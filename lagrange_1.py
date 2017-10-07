import numpy as np

def lagrange_1(x,y,u):
    """拉格朗日插值程序

    :x: TODO
    :y: TODO
    :u: TODO
    :returns: TODO

    """
    n=len(x)
    m = len(u)
    v = np.zeros(m)

    for i in range(n):
        w = np.ones(m)
        for j in range(n):
            if j!= i:
                w = (u)
