import sys
from scipy.interpolate import lagrange
import matplotlib.pyplot as plt
import numpy as np
from sympy import *


def adam_bashford_4(func, interv, h):

    x = Symbol('x') #Inicializa "x" como símbolo
    y = Symbol('y') #Inicializa "y" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy 
    fx = lambdify([x, y], f, modules=['numpy']) #Se inicializa la función f(x,y) 

    listaX = [2, 2.2, 2.4] #Valores iniciales x
    listaY = [1, 1.191, 1.365] #Valores inicial y

    x_temp = 2.6 #Valor temporal xk
    y_temp = 1.528 #Valor temporal yk

    cont = 3 #Contador para yk

    while(x_temp <= interv[1]): #Mientras se encuentre en el dominio [2,4] 


        listaX.append(float(x_temp)) #Añade a lista xk
        listaY.append(float(y_temp)) #Añade a lista yk

        y_temp = yk(listaX, listaY, h, cont, fx) #Llama función que calcula yk
        x_temp += h #Incrementa valor de x en h

        cont += 1 #Incrementa contador

    p = lagrange(listaX,listaY) #Crea polinomio de interpolación
    print("Polinomio de interpolación: \n", p) #Imprime polinomio de interpolación
    x = np.linspace(2,4,40) #Crea una lista con n valores entre el dominio [2,4]

    #Formato de la línea y los puntos
    plt.scatter(x, p(x),color='r',zorder=1) 
    plt.plot(x, p(x),color='b',zorder=2)

    #Especifiaciones de gráfica
    plt.title("Gráfica polinomio de interpolación")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.ylim([-5,5])

    #Muestra gráfica
    plt.show()
    
    return [listaX, listaY]

def yk(listaX, listaY, h, cont, fx): 

    #Función que crea yk
    yk = listaY[cont] + (h/24) * (55 * fx(listaX[cont], listaY[cont]) - 59 * fx(listaX[cont-1], listaY[cont-1]) + 37 * fx(listaX[cont-2], listaY[cont-2]) - 9 * fx(listaX[cont-3], listaY[cont-3]))

    return yk

adam_bashford_4('1 + (x - y) ** 2', [2,4], 0.2)
