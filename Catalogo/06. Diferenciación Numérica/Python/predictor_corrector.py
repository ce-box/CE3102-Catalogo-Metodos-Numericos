import sys
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
import numpy as np
from sympy import *

def predictor_corrector(func, y0, interv, N_puntos):

    x  = Symbol('x') #Inicializa "x" como símbolo
    y  = Symbol('y') #Inicializa "y" como símbolo
    f  = sympify(func) #Se traduce el string "func" a una función de sympy 
    fxy = lambdify([x, y], f, modules=['numpy']) #Se inicializa la función f(x,y) 
    
    a = interv[0]
    b = interv[1]

    listaX = [] 
    listaY = []

    h = (b-a)/(N_puntos-1)
    xn = 0
    yn = y0
    n = 0

    listaX.append(xn)
    listaY.append(yn)

    while (xn < b):
        yn_euler = euler(fxy, yn, xn, h)

        n += 1
        xn_sig = a + n * h

        yn = yn + h * (fxy(xn, yn) + fxy(xn_sig, yn_euler))/2
        xn = xn_sig

        listaX.append(round(xn, 2))
        listaY.append(round(yn, 4))

    p = lagrange(listaX,listaY) #Crea polinomio de interpolación
    print("Polinomio de interpolación: \n", p) #Imprime polinomio de interpolación
    x = np.linspace(a,b,N_puntos) #Crea una lista con n valores entre el dominio [2,4]

    #Formato de la línea y los puntos
    plt.scatter(x, p(x),color='r',zorder=1) 
    plt.plot(x, p(x),color='b',zorder=2)

    #Especifiaciones de gráfica
    plt.title("Gráfica polinomio de interpolación")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.ylim([y0-1,yn+1])

    #Muestra gráfica
    plt.show()
    
    return [listaX, listaY]

def euler(fxy, yn, xn, h):

    yn_euler = yn + h * fxy(xn,yn)

    return yn_euler 


#Valores iniciales 
if __name__ == '__main__':
    func = 'y - x^2 + 1' #Funcion
    interv = [0,2] #Intervalo 
    y0 = 0.5       #Valor inicial en y
    N_puntos = 11  #Cantidad de puntos
    predictor_corrector(func,y0,interv,N_puntos)
