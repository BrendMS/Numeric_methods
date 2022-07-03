from re import A
import numpy as np
import math
import ep2
import ep1

################ EP3 ##################

def phi(x, x0, x1, h):
    if x >= x0 & x <= x1:
        return (x-x0)/h
    else:
        return 0

def phil(x, x0, x1, h):
    if x >= x0 & x <= x1:
        return (x-x0)/h
    else:
        return 0
    
def montarMatrizA_k1_q0(n, h, x): #MATRIZ QUANDO K(X) = 1 e Q(X) = 0
    A = np.zeros(n,n)
    for i in range(0, n):
       for j in range(0, n): 
            #diagonal principal
            if i == j:
                A[i][j] = 2/h
            #diagonal superior
            elif (j - i) == 1:
                A[i][j] = -1/h
            #diagonal inferior
            elif (j - i) == -1:
                A[i][j] = -1/h
    return A

def produtointerno_f_phi(k, q, u, v, ul, vl):
    pass

def montarMatrizB(n): #VETOR SOLUCAO DA MATRIZ A
    B = np.zeros(n)
    for i in range(0,n):
        B[i] = produtointerno_f_phi()
    return B




################ FUNCOES ##################
def fx_validacao(x):
    return 12*x*(1 - x) - 2

################ MAINS ###################
def main():
    n = 3
    h = 1/(n+1)
    x = np.zeros(n+1)
    print("Desisto!")

def main_validacao():
    n = 7
    h = 1/(n+1)
    x = np.zeros(n+1)
    for i in range(0, n+1):
        x[i] = i*h
    print(x)
    #montarMatrizA_k1_q0(n, h, x)

if __name__ == "__main__":
    main_validacao()

    