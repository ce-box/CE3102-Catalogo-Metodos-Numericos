import numpy as np
import sys
from sympy import *
from matplotlib import pyplot as plt

def punto_fijo(func, rango, tol, iterMax):

    #Funcion de prueba: "E**(-x)" ó "log(3*x)", [1,2], 10**(-10), 100

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy
    df = f.diff(x) #Se calcula la derivada de "f"
    
    gx = lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)

    a = rango[0]
    b = rango[1]

    D = []
    A = []
    
    cont = 1 #Se inicializa contador

    try:
        puntosC = solve(df,x)
    except:
        print("No se garantizan puntos fijos para la función " + func)

    if(a <= gx(a) <= b and a <= gx(b) <= b):
        if(len(puntosC) == 0):
            print("No se garantizan puntos fijos para la función " + func)

        elif(len(puntosC) == 1):
            if(a <= gx(puntosC[0]) <= b):
                print("La función " + func + " tiene al menos un punto fijo")

        else:
            if(a <= gx(puntosC[0]) <= b and a <= gx(puntosC[1]) <= b):
                print("La función " + func + " tiene al menos un punto fijo")
    else:
        print("No se garantizan puntos fijos para la función " + func)

    b = gx(a)
    tramo = abs(b-a)
    
    while(cont <= iterMax):

        a = b
        b = gx(a)
        tramo = abs(b-a)

        D.append(cont)
        A.append(tramo)
        
        if(tramo <= tol):
            break
    
        cont += 1

    plt.plot(np.array(D), np.array(A))
    plt.title("Gráfico")
    plt.legend(["Error"])
    plt.show()

    return [b, tramo]
    

    
