from re import A
import numpy as np
import matplotlib.pyplot as plt
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
    elif n == 3: ##EQUILIBRIO FORÇANTES DE CALOR 4.3 - COM CONSTANTES
        if b == 1: 
            return  0#u(x)
        else:
            Qmais =  37.5  * 1000000
            Qmenos = 11.25 * 1000000
            return  Qmais - Qmenos#f(x) 
    elif n == 4: ##EQUILIBRIO FORÇANTES DE CALOR 4.3 - NAO CONSTANTE
        if b == 1: 
            return  0#u(x)
        else:
            Qmais = 7.5 * 1000000
            Qmenos = 11.25 * 1000000
            L = 0.03
            sigma = 1
            theta = 0.01
            return  (Qmais)*math.e**( (-(x-(L/2))**2) / (sigma**2) ) - (Qmenos)*( math.e**((-x**2)/(theta**2)) + math.e**((-(x-L)**2)/(theta**2)) )#f(x) 

################ MAINS ###################

def main_doismateriais(n, plotar, ks, ka, d):
    L =  0.02
    va = 20 + 273.5
    vb =  20 + 273.5
    def k(x):
        if (L/2 - d) <= x <= (L/2 + d):
            return ks
        else:
            return ka
    def q(x):
        return 0
    funcaoselect = 4
    h = L/(n+1)
    x = np.zeros(n+2)
    for i in range(0, n+2, 1):
        x[i] = (i)*h
    A, a, b, c = montarMatrizA_k1_q0(n, h, x, k, q)
    B = montarMatrizB(n, x, funcaoselect, h)
    l, u = ep1.decompLU(a ,b ,c, n)
    alphas = ep1.solucaoLU(l, u, c, B, n)

    def ubarra(xvar, alphas, x, h, n):
        ubarra = 0.0
        for j in range(1,n+1):
            ubarra += alphas[j-1]*phi(xvar,x[j-1],x[j+1],x[j],h)
        return ubarra

    resultado = np.zeros(n+2)
    for i in range(0,n+2):
        resultado[i] = ubarra(x[i], alphas, x, h, n) + va + ((vb-va) * x[i]) - 273.5 #ubarra(x[i], alphasGER, x, h, n) + va + (vb-va) * x[i] - ubarra(x[i], alphasDIS, x, h, n) - 273.5

    # RESULTADOS:
    for i in range(0,n+2):
        print("A temperaatura em x = " + str(x[i]) + "mm vale: " +  str(resultado[i]) + " °C")
    if plotar == 1:
        plt.plot(x*1000, resultado)
        plt.title('Temperatura')
        plt.ylabel('Temperatura (C°)')
        plt.xlabel('x(mm)')
        plt.show()

def main_equilibriocomforcantesdecalor(n, plotar): #4.3
    L = 0.02
    va = 20 + 273.5
    vb = 20 + 273.5
    def k(x):
        return 3.6
    def q(x):
        return 0
    #CALOR TOTAL CONSTANTE:
    funcaoselect = 4
    h = L/(n+1)
    x = np.zeros(n+2)
    for i in range(0, n+2, 1):
        x[i] = (i)*h
    A, a, b, c = montarMatrizA_k1_q0(n, h, x, k, q)
    B = montarMatrizB(n, x, funcaoselect, h)
    l, u = ep1.decompLU(a ,b ,c, n)
    alphas = ep1.solucaoLU(l, u, c, B, n)

    def ubarra(xvar, alphas, x, h, n):
        ubarra = 0.0
        for j in range(1,n+1):
            ubarra += alphas[j-1]*phi(xvar,x[j-1],x[j+1],x[j],h)
        return ubarra

    resultado = np.zeros(n+2)
    for i in range(0,n+2):
        resultado[i] = ubarra(x[i], alphas, x, h, n) + va + ((vb-va) * x[i]) - 273.5 #ubarra(x[i], alphasGER, x, h, n) + va + (vb-va) * x[i] - ubarra(x[i], alphasDIS, x, h, n) - 273.5

    # RESULTADOS:
    for i in range(0,n+2):
        print("A temperaatura em x = " + str(x[i]) + "mm vale: " +  str(resultado[i]) + " °C")
    if plotar == 1:
        plt.plot(x*1000, resultado)
        plt.title('Temperatura')
        plt.ylabel('Temperatura (C°)')
        plt.xlabel('x(mm)')
        plt.show()

def main_equilibriocomforcantesdecalor_constante(n, plotar): #4.3
    L = 0.02
    va = 20 + 273.5
    vb = 20 + 273.5
    def k(x):
        return 3.6
    def q(x):
        return 0
    #CALOR TOTAL CONSTANTE:
    funcaoselect = 3
    h = L/(n+1)
    x = np.zeros(n+2)
    for i in range(0, n+2, 1):
        x[i] = (i)*h
    A, a, b, c = montarMatrizA_k1_q0(n, h, x, k, q)
    B = montarMatrizB(n, x, funcaoselect, h)
    l, u = ep1.decompLU(a ,b ,c, n)
    alphas = ep1.solucaoLU(l, u, c, B, n)

    def ubarra(xvar, alphas, x, h, n):
        ubarra = 0.0
        for j in range(1,n+1):
            ubarra += alphas[j-1]*phi(xvar,x[j-1],x[j+1],x[j],h)
        return ubarra

    resultado = np.zeros(n+2)
    for i in range(0,n+2):
        resultado[i] = ubarra(x[i], alphas, x, h, n) + va + ((vb-va) * x[i]) - 273.5 #ubarra(x[i], alphasGER, x, h, n) + va + (vb-va) * x[i] - ubarra(x[i], alphasDIS, x, h, n) - 273.5

    #RESULTADOS:
    for i in range(0,n+2):
        print("A temperatura em x = " + str(x[i]) + "mm vale: " +  str(resultado[i]) + " °C")
    if plotar == 1:
        plt.plot(x*1000, resultado)
        plt.title('Temperatura com geração de calor constante')
        plt.ylabel('Temperatura (C°)')
        plt.xlabel('x(mm)')
        plt.show()

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
    print("EP3")
    print("Insira o valor de n:")
    n = int(input())
    print("Voce quer que os resultados sejam plotados? 1-Sim, 0-Não")
    plotar = int(input())
    print("O que vc deseja rodar?")
    print("1 - Validacao")
    print("2 - Validacao(Complemento)")
    print("3 - Equilibrio com forcantes de calor com aquecimento e resfriamento constantes")
    print("4 - Equilibrio com forcantes de calor com aquecimento e resfriamento com modelo realista")
    print("5 - Equilibrio com forcantes de calor com dois materiais no chip")
    select = int(input())
    if  select == 1:
        main_validacao(n)
    elif select == 2:
        main_validacao_comp(n)
    elif select == 3:
        main_equilibriocomforcantesdecalor_constante(n, plotar)
    elif select == 4:
        main_equilibriocomforcantesdecalor(n, plotar)
    elif select == 5:
        print("Insira o raio do chip em milimetros: (Obs: Largura total = 20mm)")
        d = float(input())/1000
        print("Insira a condutividade do chip(Silicio: 3.6): ")
        ks = float(input())
        print("Insira a condutividade do outro material: ")
        ka = float(input())
        main_doismateriais(n, plotar, ks, ka, d)
    #main_validacao(n)
    #main_validacao_comp(n)
    #main_equilibriocomforcantesdecalor_constante(n, plotar)
    #main_equilibriocomforcantesdecalor(n, plotar)
    