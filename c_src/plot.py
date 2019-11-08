import numpy as np
import matplotlib.pyplot as plt
import sys

datafiles = sys.argv[1:]

for datafile in datafiles:
    data = np.loadtxt(datafile)
    nplot, mplot = data.shape
    if nplot <= 1000:
        for i in range(1, mplot):
            plt.plot(data[:, 0], data[:, i], ".", label="plot%i" % i)
    if nplot >= 1000:
        for i in range(1, mplot):
            plt.plot(data[-1000:, 0], data[-1000:, i], ".", label="plot%i" % i)

plt.legend()
plt.show()
