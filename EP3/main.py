import numpy as np
import math

n = 3
h = 1/(n+1)
x = np.zeros(n+1)
for i in range(0, n+1):
    x[i] = i*h

def phi(x, x0, x1):
    if x >= x0 & x <= x1:
        return (x-x0)/h
    else:
        return 0

def montarMatrizA():
    pass

def main():
    print("Desisto!")

if __name__ == "__main__":
    main()

    