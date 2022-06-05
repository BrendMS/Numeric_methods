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

def transp_nos(xi, c, d):
    #CADA VALOR TRANSPORTADO
    yj = (((d-c)/2)*xi)+((c+d)/2)

    return yj

def transp_pesos(wi, c, d):
    #CADA VALOR TRANSPORTADO
    wj = ((d-c)/2)*wi

    return wj

def d(x, sel):
    #DEFINICAO DA FUNCAO DO EXTREMO SUPERIOR
    if   sel == 1:
        pass
    elif sel == 2:
        result = 1 - x**2
    elif sel == 3:
        result = x**3
    elif sel == 4:
        result = math.e**(-x**2)

    return result

def c(x, sel):
    #DEFINICAO DA FUNCAO DO EXTREMO INFERIOR
    if   sel == 1:
        #se esse for 0 também apaga essa funcao e só bota zero lá
        pass
    elif sel == 2:
        result = 0
    elif sel == 3:
        result = 0
    elif sel == 4:
        result = 0

    return result

def integraldupla(n, sel, a, b):
    #CRIA MATRIZES NECESSÁRIAS 
    xj = np.zeros(n)
    Fzao = np.zeros(n)
    y = np.zeros(n)
    v = np.zeros(n)

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

    for j in range(n):
        xj[j] = transp_nos(x[j], a, b)

    #CALCULAR F
    for i in range(0,n):
        for j in range(0,n):
            y[i][j] = transp_nos(x[i], c(xj[i], sel), d(xj[i], sel))
            v[i][j] = transp_pesos(w[i], c(xj[i], sel), d(xj[i], sel))
            Fzao[i] = Fzao[i] + v[i][j]*f(x[i],y[i][j])   

    #CALCULAR I
    I = 0
    for i in range(0,n):
        I = I + w[i]*Fzao(x[i]) 

    I = I * ((b-a)/2)
    return 

def f(x):
    return x**2

def main():
    print("Escolha qual dos exemplos do enunciado que deseja executar('1', '2', '3', '4')")
    sel = int(input())
    print("Escolha o valor do n (6, 8 ou 10)")
    n  = int(input())
    if   sel == 1:
        pass
    elif sel == 2:
        pass
    elif sel == 3:
        pass 
    elif sel == 4:
        pass 

    print("n: " + str(n))
    print("resultado da integral: " + str(result))


if __name__ == "__main__":
    main()