from re import A
import numpy as np
import math
import ep2
import ep1
################ EP2 ##################

def integralSimples(a, b, funcaoselect, x, i, h):
    resultado = 0
    w = 1 / (np.sqrt(3))
    for i in range(1,2,1):
        if i==1:
            resultado = resultado + fvezesphi( (a + ((b-a)/2) * (1-w)) , x, i, funcaoselect, h)
        else:
            resultado = resultado + fvezesphi( (a + ((b-a)/2) * (1+w)) , x, i, funcaoselect, h)
    return (resultado*((b-a)/2))


################ EP3 ##################

def phi(x, x0, x1, h):
    if (x >= x0):
        if (x <= x1):
            return ((x-x0)/h)*1.0
        else:
            return 0.0
    else:
        return 0.0

def fvezesphi(xvar, x, i, funcaoselect, h):
    return (funcaoescolhida(xvar, funcaoselect) * phi(xvar, x[i-1], x[i+1], h))
    
def montarMatrizA_k1_q0(n, h): #MATRIZ QUANDO K(X) = 1 e Q(X) = 0
    A = np.zeros((n,n))
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

def produtointerno_f_phi(x, i, funcaoselect, h):
    return integralSimples(x[i-1], x[i+1], funcaoselect, x, i, h) #integralSimples(a, b, funcaoselect, x, i):
    

def montarMatrizB(n, x, funcaoselect, h): #VETOR SOLUCAO DA MATRIZ A
    B = np.zeros(n)
    for i in range(0,n):
        B[i] = produtointerno_f_phi(x, i, funcaoselect, h)
    return B


################ FUNCOES ##################
def funcaoescolhida(x, n):
    if   n == 1:
        return 12*x*(1 - x) - 2

################ MAINS ###################
def main():
    n = 3
    h = 1/(n+1)
    x = np.zeros(n+1)
    print("Desisto!")

def main_validacao():
    funcaoselect = 1
    n = 7
    h = 1/(n+1)
    x = np.zeros(n+1)
    for i in range(0, n+1):
        x[i] = i*h
    print("VETOR X:")
    print(x)
    A = montarMatrizA_k1_q0(n, h)
    print("MATRIZ A:")
    print(A)
    B = montarMatrizB(n, x, funcaoselect, h)
    print("MATRIZ B:")
    print(B)

if __name__ == "__main__":
    main_validacao()

    