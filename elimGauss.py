import numpy as np
import sys
import time

def elim_gauss():

    a = np.array([[2, 1, 0, 0, 0, 1, 1], [1, 2, 1, 0, 0, 0, 1], [0, 1, 2, 1, 0, 0, 1], [0, 0, 1, 2, 1, 0, 1], [0, 0, 0, 1, 2, 1, 1], [1, 0, 0, 0, 1, 2, 1]])
    n = len(a)
    A = np.copy(a[:,:-1])
    b = np.copy(a[:,n]).reshape((n))
    
    #Guardo tanto L quanto U em uma unica matriz!!!!
    LU = np.eye(n)
    # @ faz multiplicacao de matrizes usando numpy
    for i in range(n):
        #Varre linhas superiores (Upper)
        LU[i,i:] = A[i,i:]-LU[i,:i] @ LU[:i,i:]
        #Varre colunas inferiores (Lower)
        LU[(i+1):,i] = ( A[(i+1):,i]- LU[(i+1):,:i] @ LU[:i,i] ) / LU[i,i]
    
    print("LU: ", LU)
    # Substituicao
    # LUx=b
    # Ly=b
    y = np.zeros(n)
    # Ly=b
    y[0] = b[0]
    for i in range(1,n,1):
        #Vetorizei aqui!
        y[i] = (b[i] - np.dot(LU[i,:i], y[:i]))
    
    # Ux=y
    x = np.zeros(n)
    x[n-1] = y[n-1]/LU[n-1,n-1]
    for i in range(n-2,-1,-1):
        #Vetorizei aqui!
        x[i] = (y[i] - np.dot(LU[i,i+1:], x[i+1:]))/LU[i,i]
    
    #Testes
    U = np.triu(LU)
    L = np.tril(LU)
    np.fill_diagonal(L, 1.0)
    print("L: ", L)
    print("U: ", U)

    return x, y

if __name__ == "__main__":
    elim_gauss()
