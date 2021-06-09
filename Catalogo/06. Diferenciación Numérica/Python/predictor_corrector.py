import sys
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
import numpy as np
from sympy import *

#Parametros de entrada
#func: funcion de dos variables, y0: valor inicial, interv: interbalo del funcion
# N_puntos: cantidad de puntos 
#Salida
#resultado: polinomio de interpolacion
#Funcion Predictor corrector
def predictor_corrector(func, y0, interv, N_puntos):

    x  = Symbol('x') #Inicializa "x" como símbolo
    y  = Symbol('y') #Inicializa "y" como símbolo
    f  = sympify(func) #Se traduce el string "func" a una función de sympy 
    fxy = lambdify([x, y], f, modules=['numpy']) #Se inicializa la función f(x,y) 
    
    #Variables a  y b
    a = interv[0]
    b = interv[1]

    #Listas de puntos
    listaX = [] 
    listaY = []

    #Valor de h, X0, Y0, n
    h = (b-a)/(N_puntos-1)
    xn = 0
    yn = y0
    n = 0

    listaX.append(xn)
    listaY.append(yn)

    #Clico principal de ejecucion 
    while (xn < b):
        #Calculo de yn+1 con Euler
        yn_euler = euler(fxy, yn, xn, h)

        n += 1 #n siguiente
        xn_sig = a + n * h #xn+1

        #Calculo de yn+1 predictor corrector
        yn = yn + h * (fxy(xn, yn) + fxy(xn_sig, yn_euler))/2
        xn = xn_sig

        listaX.append(round(xn, 2)) #Punto a la lista
        listaY.append(round(yn, 4)) #Punto a la lista

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
    
    #Resultado
    return p

#Parametros de entrada
#fxy: funcion de dos variables, yn: valor inicial, xn: valor inicial, h: valor h
#Salida
#resultado: yn+1 por euler
#Funcion Euler
def euler(fxy, yn, xn, h):
    
    #Calculo yn+1 con euler
    yn_euler = yn + h * fxy(xn,yn)

    #Resultado
    return yn_euler 


#Valores iniciales 
if __name__ == '__main__':
    func = 'y - x^2 + 1' #Funcion
    interv = [0,2] #Intervalo 
    y0 = 0.5       #Valor inicial en y
    N_puntos = 11  #Cantidad de puntos
    predictor_corrector(func,y0,interv,N_puntos)
