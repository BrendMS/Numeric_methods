import numpy as np
import matplotlib.pyplot as plt

'''
Implemente o algoritmo descrito acima para a decomposicao LU de uma matriz tridiagonal A n × n.
As matrizes A, L e U devem ser armazenadas em vetores conforme descrito no texto. Implemente
tambem o algoritmo para a resolucao de um sistema linear tridiagonal usando a decomposicao LU da
matriz. Faca as implementacoes de forma que elas possam ser usadas como partes de outros programas.
'''

def montaMatriz(n): #Função responsável por gerar os vetores 
    #a,b,c,d,bT,cT,v,w
    d  = np.zeros(n)
    a  = np.zeros(n)
    b  = np.zeros(n)
    c  = np.zeros(n)
    bT = np.zeros(n)
    cT = np.zeros(n)
    aT  = np.zeros(n)
    dT = np.zeros(n)    
    v  = np.zeros(n-1)

    for i in range(n):
        if i == (n-1): #Gerar vetor a
            a[i] = ((2*n)-1)/(2*n)
        else:
            a[i] = ((2*(i+1))-1)/(4*(i+1))

        b[i] = 2 #Gerar vetor b

        c[i] = 1 - a[i] #Gerar vetor c
        d[i] = np.cos((2*np.pi*(i+1)*(i+1))/(n*n)) #Gerar vetor d

    #Gerar vetores da matriz tridiagonal 
    bT = np.copy(b[:n-1])
    aT = np.copy(a[:n-1])
    cT = np.copy(c[:n-1])

    aT[0] = 0
    cT[-1] = 0 

    dT = np.copy(d[:n-1]) #Gerar vetor d_tio  

    v[0] = a[0] #Gerar vetor v
    v[-1] = c[n-2] 

    print(aT)
    print(bT)
    print(cT)
    print(dT)
    print(a)
    print(b)
    print(c)
    print(d)
    print(v)
    
    return a,b,c,d,aT, bT, cT, dT, v #Retornar todos os vetores dessa função

def decompLU(a, b, c, n): #Função responsável por achar os valores das matrizes L e U
    l = np.zeros(n)
    u = np.zeros(n)

    u[0] = b[0]
    for i in range(1,n) :
        if not (u[i-1] == 0):
            l[i] = a[i]/u[i-1] #Achar valores de L
        else:
            l[i] = 0

        u[i] = b[i] - l[i]*c[i-1] #Achar valores de U
        

    #print("A matriz L: " + str(l))
    #print("A matriz U: " + str(u))
    return l, u
            
def solucaoLU(l, u, c, d, n):
    y = np.zeros(n)

    #### SOLUCAO DO SISTEMA DE L
    y[0] = d[0]
    for i in range(1, n):
        y[i] = d[i]-l[i]*y[i-1]
    #print("A solucao do sistema L(y) obtida: " + str(y))

    #### SOLUCAO DO SISTEMA DE U
    x = np.zeros(n)
    x[n-1] = y[n-1]/u[n-1] 
    for i in range(n-2, -1, -1):
        if not (u[i] == 0):
            x[i] = (y[i]-c[i]*x[i+1])/u[i]
        else:
            x[i] = 0
    #print("A solucao do sistema U(x) obtida: " + str(x))
    return x

def acharxT(zT, yT, a, b, c, d, n):
    xN = (d[n-1] - c[n-1]*yT[0] - a[n-1]*yT[n-2]) / (b[n-1] - c[n-1]*zT[0] -a[n-1]*zT[n-2])
    xT = yT - xN*zT
    xT = np.append(xT, xN)

    return xT
    
def resolveciclica(n):
    a, b, c, d, aT, bT, cT, dT, v = montaMatriz(n)

    l_c, u_c = decompLU(aT, bT, cT, n-1)
    y = solucaoLU(l_c, u_c, cT, dT, n-1)
    z = solucaoLU(l_c, u_c, cT, v, n-1)

    result = acharxT(z, y, a, b, c, d, n)

    return result

def main():
    print("Qual o tamanho da matriz nxn? (Mínimo 3x3)") #FAZER INPUT DA MATRIZ E DOS 4 VETORES TAMBEM
    n = int(input())
    #print("Valor selecionado para n: " + str(n))

    xT = resolveciclica(n)
    print("O resultado é:")
    print(xT)

if __name__ == "__main__":
    main()