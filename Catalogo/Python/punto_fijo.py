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

    D = [] #Lista para graficar eje x (iteraciones)
    A = [] #Lista para graficar eje y (errores)
    
    cont = 1 #Se inicializa contador

    #Verifica si se puede garantizar un punto fijo para la función ingresada

    try:
        puntosC = solve(df,x)
    except:
        print("No se garantizan puntos fijos para la función " + func)

    if(a <= gx(a) <= b and a <= gx(b) <= b):
        if(len(puntosC) == 0): #Si no hay puntos críticos
            print("No se garantizan puntos fijos para la función " + func)

        elif(len(puntosC) == 1): #Si solo hay un punto crítico
            if(a <= gx(puntosC[0]) <= b):
                print("La función " + func + " tiene al menos un punto fijo")

        else: #Si hay dos puntos críticos
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

        D.append(cont) #Se añade a la lista de iteraciones
        A.append(tramo) #Lista de errores
        
        if(tramo <= tol):
            break
    
        cont += 1 #Incrementa contador

    plt.plot(np.array(D), np.array(A))
    plt.title("Gráfico")
    plt.legend(["Error"])
    plt.show()

    return [b, tramo] #Retorna la aproximación y el error
    

    
