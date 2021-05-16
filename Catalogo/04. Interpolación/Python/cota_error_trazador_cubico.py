import numpy as np
import sys
from sympy import *

import trazador_cubico as tc

def cota_traz_cubico(func, numPuntos):

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy 
    fx = lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)

    d1f = f.diff(x)
    d2f = d1f.diff(x)
    d3f = d2f.diff(x)
    d4f = d3f.diff(x)
    d4fx = lambdify(x, d4f, modules=['numpy']) #Se inicializa la función f(x)
    

    Xk = [1,1.05,1.07,1.1]
    Yk = [2.718282, 3.286299, 3.527609, 3.905416]
    
    trazCub = tc.traz_cubico(Xk, Yk)

    functCont = 0
    h = max(trazCub[1])
    
    while(functCont < len(trazCub[0])):
        
        numCont = 0
        
        while(numCont < numPuntos):
            
            sx = lambdify(x, trazCub[0][functCont], modules=['numpy'])
            
            error = abs(fx(numCont) - sx(numCont))
            comp = (5 * h**4)*d4fx(numCont)/384

            print("\n\n" + str(error))
            print(str(comp) + "\n\n")
            
            if(error < comp):
                print("La cota del error es:" + str(error))
                return comp

            numCont += 1
        
        
        functCont += 1

cota_traz_cubico('3*x*E**x - 2*E**x', 10)
