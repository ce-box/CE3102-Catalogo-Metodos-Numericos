import numpy as np
import sys
from sympy import *
from matplotlib import pyplot as plt

def newton_raphson(func, xk, tol, iterMax):

    #Funcion de prueba: "cos(2*x)**2 - x**2", 3/4, 10**(-10), 100

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy
    df = f.diff(x) #Se calcula la derivada de "f"
    
    fx = lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)
    dfx = lambdify(x, df, modules=['numpy']) #Se inicializa la derivada de la función f(x)

    D = []
    A = []

    cont = 1
    #Se inicializa el contador

    #Verifica que la tolerancia no sea negativa
    if (tol < 0):
        return print("La tolerancia debe ser mayor que cero")

    while(cont < iterMax):

        if(dfx(xk) == 0): #Indefinición de denominador
            break

        xk = xk - ((fx(xk))/(dfx(xk)))
        y = fx(xk)

        D.append(cont)
        A.append(y)

        if(np.absolute(y) < tol):
            break

        cont += 1 #Incrementa contador

    plt.plot(np.array(D), np.array(A))
    plt.title("Gráfico")
    plt.legend(["Error"])
    plt.show() 

    return [xk, np.absolute(y)] #Retorna la aproximación del cero de la función, el error, y la cantidad de iteraciones
    
