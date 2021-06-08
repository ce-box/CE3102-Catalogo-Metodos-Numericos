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

    listaX = [2, 2.2, 2.4]
    listaY = [1, 1.191, 1.365]

    x_temp = 2.6
    y_temp = 1.528

    cont = 3

    while(x_temp <= interv[1]):      


        listaX.append(float(x_temp))
        listaY.append(float(y_temp))

        y_temp = yk(listaX, listaY, h, cont, fx)
        x_temp += h

        cont += 1

    p = lagrange(listaX,listaY)
    print("Polinomio de interpolación: ", p)
    x = np.linspace(2,4,40)
    plt.scatter(x, p(x),color='r',zorder=1)
    plt.plot(x, p(x),color='b',zorder=2)

    plt.title("Gráfica polinomio de interpolación")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.ylim([-5,5])
    
    plt.show()
    
    return [listaX, listaY]

def yk(listaX, listaY, h, cont, fx):

    yk = listaY[cont] + (h/24) * (55 * fx(listaX[cont], listaY[cont]) - 59 * fx(listaX[cont-1], listaY[cont-1]) + 37 * fx(listaX[cont-2], listaY[cont-2]) - 9 * fx(listaX[cont-3], listaY[cont-3]))

    return yk

resultado = adam_bashford_4('1 + (x - y) ** 2', [2,4], 0.2)
