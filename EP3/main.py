from re import A
import numpy as np
import math
import ep2
import ep1
################ EP2 ##################

def integralSimples_fphi(a, b, funcaoselect, x, i, h):
    w = np.sqrt(3)/3
    t = (w/2)*(b-a)
    t1 = (a+b)/2
    resultado = fvezesphi(t + t1, x, i, funcaoselect, h) + fvezesphi(-t + t1, x, i, funcaoselect, h) 
    return resultado*((b-a)/2)

def integralSimples_phiphi(a, b, x, j, i, h, k, q):
    w = np.sqrt(3)/3
    t = (w/2)*(b-a)
    t1 = (a+b)/2
    resultado = phiphi(t + t1, k, q, x, i, j, h) + phiphi(-t + t1, k, q, x, i, j, h) #(xvar, k, q, x, i, j, h)
    return resultado*((b-a)/2)


################ EP3 ##################

def phi(x, x0, x1, xi, h):
    if  (x0 <= x <= xi):
        return ((x-x0)/h)
    elif (xi <= x <= x1):
        return ((x1-x)/h)
    else:
        return 0

def phi_l(x, x0, x1, xi, h):
    if  (x0 <= x <= xi):
        return (1/h)
    elif (xi <= x <= x1):
        return (-1/h)
    else:
        return 0


def fvezesphi(xvar, x, i, funcaoselect, h):
    return (funcaoescolhida(xvar, funcaoselect, 0) * phi(xvar, x[i-1], x[i+1], x[i], h))

def phiphi(xvar, k, q, x, i, j, h):
    phiphi = ( k(xvar)*phi_l(xvar, x[i-1], x[i+1], x[i], h)*phi_l(xvar, x[j-1], x[j+1], x[j], h) + q(xvar)*phi(xvar, x[i-1], x[i+1], x[i], h)*phi(xvar, x[j-1], x[j+1], x[j], h) )
    return phiphi

def produtointerno_phiphi(q, k, x, j, i, h):
    return integralSimples_phiphi(x[i-1], x[i], x, j, i, h, k, q) + integralSimples_phiphi(x[i], x[i+1], x, j, i, h, k, q)#(L, x, j, i, h, k, q):

def produtointerno_f_phi(x, i, funcaoselect, h):
    return ((integralSimples_fphi(x[i-1], x[i], funcaoselect, x, i, h) + integralSimples_fphi(x[i], x[i+1], funcaoselect, x, i, h))) #integralSimples(a, b, funcaoselect, x, i):
    
def montarMatrizA_k1_q0(n, h, x, k, q): #MATRIZ QUANDO K(X) = 1 e Q(X) = 0
    A = np.zeros((n,n))
    Am = np.zeros(n)
    As = np.zeros(n)
    Ai = np.zeros(n)
    for i in range(0, n):
       for j in range(0, n): 
            #diagonal principal
            if i == j: #2/h
                valor = produtointerno_phiphi(q, k, x, j+1, i+1, h)
                A[i][j] = valor
                Am[j] = valor
            #diagonal superior
            elif (j - i) == 1: #-1/h
                valor = produtointerno_phiphi(q, k, x, j+1, i+1, h)
                A[i][j] = valor
                As[j] = valor
            #diagonal inferior
            elif (j - i) == -1: #-1/h
                valor = produtointerno_phiphi(q, k, x, j+1, i+1, h)
                A[i][j] = valor
                Ai[j] = valor
    return A, As, Am, Ai 
    

def montarMatrizB(n, x, funcaoselect, h): #VETOR SOLUCAO DA MATRIZ A
    B = np.zeros(n)
    for i in range(0,n):
        B[i] = produtointerno_f_phi(x, i+1, funcaoselect, h)
    return B

def calcularErro(n, funcaoselect, alphas, h, x, it):

    def ubarra(xvar, alphas, x, h, n):
        ubarra = 0.0
        for j in range(1,n+1):
            ubarra += alphas[j-1]*phi(xvar,x[j-1],x[j+1],x[j],h)
        return ubarra

    erromax = 0.0
    xvar = np.linspace(0,1,it)

    for i in range(0,it):
        U = funcaoescolhida(xvar[i], funcaoselect, 1)
        Ubarra = ubarra(xvar[i], alphas, x, h, n)
        erro = abs(Ubarra - U)
        if erro >= erromax:
            erromax = erro

    return erromax

    

################ FUNCOES ##################
def funcaoescolhida(x, n, b):
    if   n == 1: ## VALIDACAO
        if b == 1: 
            return (x**2)*((x-1)**2) #u(x)
        else:
            return (12*x*(1 - x)) - 2 #f(x)

    elif n == 2: ##COMPLEMENTO
        if b == 1: 
            return ((x-1)*(math.e**(-x) - 1)) #u(x)
        else:
            return (math.e**x + 1) #f(x)
    elif n == 3: ##EQUILIBRIO FORÇANTES DE CALOR 4.3
        if b == 1: 
            return  0#u(x)
        else:
            return  (37500000)*math.e**( (-(x-0.02/2)**2) / (0.5**2) )#f(x) 


################ MAINS ###################

def main_equilibriocomforcantesdecalor(n): #4.3
    L = 0.02
    def k(x):
        return 3.6
    def q(x):
        return 0
    funcaoselect = 3
    h = L/(n+1)
    x = np.zeros(n+2)
    for i in range(0, n+2, 1):
        x[i] = (i)*h
    #print("VETOR X:")
    #print(x)
    A, a, b, c = montarMatrizA_k1_q0(n, h, x, k, q)
    #print("MATRIZ A (As(a), Am(b), Ai(c)):")
    #print(A)
    #print(a)
    #print(b)
    #print(c)
    B = montarMatrizB(n, x, funcaoselect, h)
    #print("MATRIZ B(d):")
    #print(B)

    ##DECOMPOSIÇÃO LU
    l, u = ep1.decompLU(a ,b ,c, n)
    alphas = ep1.solucaoLU(l, u, c, B, n)
    #print("ALPHAS / Solucao do LU:")
    #print(alphas)

    def ubarra(xvar, alphas, x, h, n):
        ubarra = 0.0
        for j in range(1,n+1):
            ubarra += alphas[j-1]*phi(xvar,x[j-1],x[j+1],x[j],h)
        return ubarra

    x0 = ubarra(0, alphas, x, h, n)
    x1 = ubarra(0.010, alphas, x, h, n)
    x2 = ubarra(0.015, alphas, x, h, n)
    x3 = ubarra(0.020, alphas, x, h, n)

    print("Valor de calor em x =  0mm: " + str(x0))
    print("Valor de calor em x = 10mm: " + str(x1))
    print("Valor de calor em x = 15mm: " + str(x2))
    print("Valor de calor em x = 20mm: " + str(x3))


    #Calcular erro
    


def main_validacao_comp(n): # 4.2 complemento
    L = 1
    def k(x):
        return math.e**x 
    def q(x):
        return 0
    funcaoselect = 2
    h = L/(n+1)
    x = np.zeros(n+2)
    for i in range(0, n+2, 1):
        x[i] = (i)*h
    print("VETOR X:")
    print(x)
    A, a, b, c = montarMatrizA_k1_q0(n, h, x, k, q)
    print("MATRIZ A (As(a), Am(b), Ai(c)):")
    #print(A)
    print(a)
    print(b)
    print(c)
    B = montarMatrizB(n, x, funcaoselect, h)
    print("MATRIZ B(d):")
    print(B)

    ##DECOMPOSIÇÃO LU
    l, u = ep1.decompLU(a ,b ,c, n)
    alphas = ep1.solucaoLU(l, u, c, B, n)
    print("ALPHAS / Solucao do LU:")
    print(alphas)

    #Calcular erro
    erro = calcularErro(n, funcaoselect, alphas, h, x, 1000)
    print("Erro encontrado: " + str(erro))

def main_validacao(n): # 4.2
    L = 1
    def k(x):
        return 1
    def q(x):
        return 0
    funcaoselect = 1
    h = L/(n+1)
    x = np.zeros(n+2)
    for i in range(0, n+2, 1):
        x[i] = (i)*h
    print("VETOR X:")
    print(x)
    A, a, b, c = montarMatrizA_k1_q0(n, h, x, k, q)
    print("MATRIZ A (As(a), Am(b), Ai(c)):")
    #print(A)
    print(a)
    print(b)
    print(c)
    B = montarMatrizB(n, x, funcaoselect, h)
    print("MATRIZ B(d):")
    print(B)

    ##DECOMPOSIÇÃO LU
    l, u = ep1.decompLU(a ,b ,c, n)
    alphas = ep1.solucaoLU(l, u, c, B, n)
    print("ALPHAS / Solucao do LU:")
    print(alphas)

    #Calcular erro
    erro = calcularErro(n, funcaoselect, alphas, h, x, 1000)
    print("Erro encontrado: " + str(erro))

if __name__ == "__main__":
    n = 63
    #main_validacao(n)
    #main_validacao_comp(n)
    main_equilibriocomforcantesdecalor(n)
    