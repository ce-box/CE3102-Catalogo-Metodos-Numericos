import numpy as np
import sys
from sympy import *

import trazador_cubico as tc

def cota_traz_cubico(func, numPuntos):

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy 
    fx = lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)

    d1f = f.diff(x) #Se calcula la primera derivada de f
    d2f = d1f.diff(x) #Se calcula la segunda derivada de f
    d3f = d2f.diff(x) #Se calcula la tercera derivada de f
    d4f = d3f.diff(x) #Se calcula la cuarta derivada de f
    d4fx = lambdify(x, d4f, modules=['numpy']) #Se inicializa la función de la cuarta derivada de f
    

    Xk = [1,1.05,1.07,1.1] #Valores iniciales de x
    Yk = [2.718282, 3.286299, 3.527609, 3.905416] #Imágenes de los valores inciales de x
    
    trazCub = tc.traz_cubico(Xk, Yk) #Se llama a la función traz_cubico para conocer S(X)

    functCont = 0 #Contador para funciones sn(x)
    h = max(trazCub[1]) #Valor mas alto de h
    
    while(functCont < len(trazCub[0])):
        
        numCont = 0 #Contador de puntos
        
        while(numCont < numPuntos):
            
            sx = lambdify(x, trazCub[0][functCont], modules=['numpy']) #Se crea la funcion s(x) en cont
            
            error = abs(fx(numCont) - sx(numCont)) #Se calcula la cota de error
            comp = (5 * h**4)*d4fx(numCont)/384 #Se calcula la comparacion de la cota

            print("\n\n" + str(error))
            print(str(comp) + "\n\n")
            
            if(error < comp): #Si la cota es menor que la el maximo error aceptado
                print("La cota de error del trazador cubico natural es:" + str(error))
                return comp
 
            numCont += 1 #Incrementa contador de puntos
        
        
        functCont += 1 #Incrementa contador de s(x)

cota_traz_cubico('3*x*E**x - 2*E**x', 10)
