#!/bin/env python
from sympy import *


def main():
    x,y,z=symbols("x:z")
    x0 = symbols("x0")
    print(integrate(4/(1+x**2),x))
    print(integrate(1/(cos(x)-cos(x0)),(x,0,x0)))


if __name__ == "__main__":
    main()
