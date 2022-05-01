import numpy as np

'''
Implemente o algoritmo descrito acima para a decomposicao LU de uma matriz tridiagonal A n × n.
As matrizes A, L e U devem ser armazenadas em vetores conforme descrito no texto. Implemente
tambem o algoritmo para a resolucao de um sistema linear tridiagonal usando a decomposicao LU da
matriz. Faca as implementacoes de forma que elas possam ser usadas como partes de outros programas.
'''

def inputMatriz(n):
    matriz = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            print("De o valor " + str(i+1) + "," + str(j+1))
            matriz[i][j] = input()
            
    print("Matriz A:")
    print(str(matriz))

    #INPUT DO VETOR RESPOSTA
    print("Agora faça o input do vetor resposta")
    d = np.zeros(n)
    for i in range(n):
        print("De o valor " + str(i + 1))
        d[i] = input()
    
    print("Vetor resposta;")
    print(str(d))

    #Separar vetores a, b e c
    a = np.zeros(n)
    b = np.zeros(n)
    c = np.zeros(n)

    a[0] = matriz[0][n-1]
    c[n-1] = matriz[n-1][0]
    for i in range(n):
        b[i] = matriz[i][i]
        if (i != (n-1)):
            a[i+1] = matriz[i+1][i]
        if (i != (n-1)):
            c[i] = matriz[i][i+1]

    print("Vetor A:")
    print(str(a))
    print("Vetor B:")
    print(str(b))
    print("Vetor C:")
    print(str(c))

    return a,b,c,d

def montaMatrizN(n): #Função responsável por gerar os vetores 
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

    v[0] = a[0] #Gerar vetor v
    v[-1] = c[n-2] 

    return a,b,c,d,v #Retornar todos os vetores dessa função

def montarArray(x): 
    y = np.zeros(len(x))
    n = len(x)

    for i in range(0, len(x)):
        if x[i] != "[" and x[i] != "]" and x[i] !=",":
            y[i] = int(x[i])
    
    return y, n

def montarMatrizVet():
    
    print("qual o vetor a? (modelo: [1,2,3]")
    a = str(input())
    aV, n = montarArray(a)

    print("qual o vetor b? (modelo: [1,2,3]")
    b = str(input())
    bV = montarArray(b) 

    print("qual o vetor c? (modelo: [1,2,3]")
    c = str(input())
    cV = montarArray(c) 

    print("qual o vetor d? (modelo: [1,2,3]")
    d = str(input())
    dV = montarArray(d)
    
    v = np.zeros(n-1)
    v[0] = a[0] #Gerar vetor v
    v[-1] = c[n-2]  

    return aV, bV, cV, dV, v, n

def tornarnaociclica(a,b,c,d,n):
    #Gerar vetores da matriz tridiagonal 
    bT = np.copy(b[:n-1])
    aT = np.copy(a[:n-1])
    cT = np.copy(c[:n-1])

    aT[0] = 0
    cT[-1] = 0 

    dT = np.copy(d[:n-1]) #Gerar vetor d_tio   

    return aT, bT, cT, dT,

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
    #monta as matrizes
    a, b, c, d, v = montaMatrizN(n)
    print("Vetor A: " + str(a))
    print("Vetor B: " + str(b))
    print("Vetor C: " + str(c))
    aT, bT, cT, dT = tornarnaociclica(a, b, c, d, n)

    #faz a decomposicao
    l_c, u_c = decompLU(aT, bT, cT, n-1)
    y = solucaoLU(l_c, u_c, cT, dT, n-1)
    z = solucaoLU(l_c, u_c, cT, v, n-1)

    #encontra o valor final
    result = acharxT(z, y, a, b, c, d, n)

    return result

def main():
    print("A matriz será cíclica? (respoder com s ou n)")
    tipo = str(input())

    print("Qual tipo de entrada?") 
    print("a) matriz")
    print("b) vetores a, b, c, d")
    print("c) n")

    entrada = str(input())

    if tipo == "s":
        if entrada == "a":
            print("Qual o tamanho da matriz A? (nxn)")
            n = int(input())
            #INPUT DA MATRIZ
            a,b,c,d = inputMatriz(n)
            v = np.zeros(n-1) #Gerar vetor v
            v[0] = a[0] 
            v[-1] = c[n-2]

            aT, bT, cT, dT = tornarnaociclica(a, b, c, d, n)

            #faz a decomposicao
            l_c, u_c = decompLU(aT, bT, cT, n-1)
            y = solucaoLU(l_c, u_c, cT, dT, n-1)
            z = solucaoLU(l_c, u_c, cT, v, n-1)

            #encontra o valor final
            xT = acharxT(z, y, a, b, c, d, n)
        
        if entrada == "b":
            a, b, c, d, v, n = montarMatrizVet()
            aT, bT, cT, dT = tornarnaociclica(a, b, c, d, n)

            #faz a decomposicao
            l_c, u_c = decompLU(aT, bT, cT, n-1)
            y = solucaoLU(l_c, u_c, cT, dT, n-1)
            z = solucaoLU(l_c, u_c, cT, v, n-1)

            #encontra o valor final
            xT = acharxT(z, y, a, b, c, d, n)
            
        if entrada == "c":
            print("Qual o tamanho da matriz A? (nxn)")
            n = int(input())
            xT = resolveciclica(n)
    
    if tipo == "n":
        if entrada == "a":
          pass  
        if entrada == "b":
            a, b, c, d, v, n = montarMatrizVet()
            
            l, u = decompLU(a, b, c, n)
            xT = solucaoLU(l, u, c, d, n)

        if entrada == "c":
            a, b, c, d, v = montaMatrizN(n)
            l, u = decompLU(a, b, c, n)
            xT = solucaoLU(l, u, c, d, n)

    print("O resultado é:")
    print(xT)


if __name__ == "__main__":
    main()