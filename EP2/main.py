import numpy as np
import math

# N=6
xn6  = np.array([0.2386191860831969086305017, 0.6612093864662645136613996, 0.9324695142031520278123016, -0.2386191860831969086305017, -0.6612093864662645136613996, -0.9324695142031520278123016])
wn6  = np.array([0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961, 0.4679139345726910473898703, 0.3607615730481386075698335, 0.1713244923791703450402961])

# N=8
xn8  = np.array([0.1834346424956498049394761, 0.5255324099163289858177390, 0.7966664774136267395915539, 0.9602898564975362316835609, -0.1834346424956498049394761, -0.5255324099163289858177390, -0.7966664774136267395915539, -0.9602898564975362316835609])
wn8  = np.array([0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314, 0.3626837833783619829651504, 0.3137066458778872873379622, 0.2223810344533744705443560, 0.1012285362903762591525314])

# N=10
xn10 = np.array([0.1488743389816312108848260, 0.4333953941292471907992659, 0.6794095682990244062343274, 0.8650633666889845107320967, 0.9739065285171717200779640, -0.1488743389816312108848260, -0.4333953941292471907992659, -0.6794095682990244062343274, -0.8650633666889845107320967, -0.9739065285171717200779640])
wn10 = np.array([0.2955242247147528701738930, 0.2692667193099963550912269, 0.2190863625159820439955349, 0.1494513491505805931457763, 0.0666713443086881375935688, 0.2955242247147528701738930, 0.2692667193099963550912269, 0.2190863625159820439955349, 0.1494513491505805931457763, 0.0666713443086881375935688])

def transp_nos(n, x, c, d):
    #TRANSPORTE DE TODOS OS VALORES DOS NÓS
    yj = np.zeros(n)
    for i in range(n):
        yj[i] = (((d-c)/2)*x[i])+((c+d)/2)

    return yj

def transp_pesos(n, w, c, d):
    #TRANSPORTE DE TODOS OS VALORES DOS PESOS
    wj = np.zeros(n)
    for i in range(n):
        wj[i] = ((d-c)/2)*w[i]

    return wj

def d(x, integr):
    #DEFINICAO DA FUNCAO DO EXTREMO SUPERIOR
    #EXEMPLO 1
    if integr == 1:
        result = 1
    elif integr == 2:
        result = 1-x

    #EXEMPLO 2
    elif integr == 3:
        result = 1 - x**2
    elif integr == 4:
        result = math.sqrt(1-x)

    #EXEMPLO 3
    elif integr == 5 or integr == 6:
        result = x**2

    #EXEMPLO 4
    elif integr == 7:
        result = (1 - (x + 0.75)**2)**0.5
    elif integr == 8:
        result = math.e**(-x**2)

    return result

def c(x, integr):
    #DEFINICAO DA FUNCAO DO EXTREMO INFERIOR
    if integr == 5 or integr == 6:
        result = x**3
    else:
        result = 0

    return result

def f(x, y, integr):
    resp = 0
    if integr == 1 or integr == 3 or integr == 4:
        resp = 1
    elif integr == 2:
        resp = 1-x-y

    elif integr == 5: 
        resp = math.sqrt( (((math.e**(y/x))*y/x**2)**2) + (((math.e**(y/x))/x)**2) + 1)
    
    elif integr == 6:
        resp = math.e**(y/x)
        
    elif integr == 7 or integr == 8:
        resp = y

    return resp

def integraldupla(n, integr, a, b):
    #CRIA MATRIZES NECESSÁRIAS 
    xj = np.zeros(n)
    wj = np.zeros(n)
    Fzao = np.zeros(n)
    w = np.zeros(n)
    x = np.zeros(n)

    s = (n, n)
    y = np.zeros(s)
    v = np.zeros(s)

    #PEGAR NÓS E PESOS
    if n == 8:
        w = wn8
        x = xn8
    elif n == 10:
        w = wn10
        x = xn10
    else:
        w = wn6
        x = xn6

    #MUDANÇAS DE VARIÁVEL EM A E B 
    xj = transp_nos(n, x, a, b)
    wj = transp_pesos(n, w, a, b)

    #CALCULAR F
    for i in range(0,n):
        y[i] = transp_nos(n, x, c(xj[i], integr), d(xj[i], integr))
        v[i] = transp_pesos(n, w, c(xj[i], integr), d(xj[i], integr))

        for j in range(0,n):
            Fzao[i] +=  v[i][j]*f(xj[i],y[i][j], integr)  

    #CALCULAR I
    I = 0
    for i in range(0,n):
        I += wj[i]*Fzao[i] 
    return I 

def main():
    print("Escolha qual dos exemplos do enunciado que deseja executar('1', '2', '3', '4')")
    sel = int(input())
    print("Escolha o valor do n (6, 8 ou 10)")
    n  = int(input())
    print("n: " + str(n))  

    if   sel == 1:
        integr = 1 #cubo
        result = integraldupla(n, integr, 0, 1)
        print("Area do cubo: " + str(result))

        integr = 2 #tetraedro
        result = integraldupla(n, integr, 0, 1)
        print("Area do tetraedro: " + str(result)) 

    elif sel == 2:
        integr = 3 # primeira integral
        result = integraldupla(n, integr, 0, 1)
        print("Area 1: " + str(result))

        integr = 4 #segunda integral
        result = integraldupla(n, integr, 0, 1)
        print("Area 2: " + str(result))  
                  
    elif sel == 3:
        integr = 5 # superfície
        result = integraldupla(n, integr, 0.1, 0.5)
        print("Área da superfície: " + str(result))

        integr = 6 #região abaixo dela
        result = integraldupla(n, integr, 0.1, 0.5)
        print("Volume da região abaixo da superfície: " + str(result))

    elif sel == 4:
        integr = 7 #Calota esférica
        result = 2*np.pi*integraldupla(n, integr, 0, 0.25)
        print("Volume da calota esférica: " + str(result))    
        
        integr = 8 #Sólido de revolução
        result = 2*np.pi*integraldupla(n, integr, -1, 1)
        print("Volume do sólido de revolução: " + str(result))    

if __name__ == "__main__":
    main()