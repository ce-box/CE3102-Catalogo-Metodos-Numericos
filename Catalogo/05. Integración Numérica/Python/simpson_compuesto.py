import sys
from sympy import *

def simpson_compuesto(func, interv):

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy 
    fx = lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)

    d1f = f.diff(x) #Se calcula la primera derivada de f
    d2f = d1f.diff(x) #Se calcula la segunda derivada de f
    d3f = d2f.diff(x) #Se calcula la tercera derivada de f
    d4f = d3f.diff(x) #Se calcula la cuarta derivada de f
    d4fx = lambdify(x, d4f, modules=['numpy']) #Se inicializa la función de la cuarta derivada de f

    a = interv[0]
    b = interv[-1]
    m = 7

    h = (b-a)/(m-1)
    x = []
    
    cont = 0

    while(cont < m):
        
        xk = a+cont*h
        x.append(xk)
        cont += 1

    cont = 1
    sumPar = 0
    sumImp = 0
            
    while(cont < (len(x) - 1)):

        if(par(cont)):
            sumPar += fx(x[cont])
        else:
            sumImp += fx(x[cont])
        cont += 1

    aprox = (h/3)*(fx(x[0]) + 2*sumPar + 4*sumImp + fx(x[-1]))
    error = (((b - a)*(h**4))/180) * abs(d4fx(2))

    print([aprox, error])
    return[aprox, error]
    

def par(x):
    if((x%2) == 0):
        return True
    else:
        return False

simpson_compuesto('log(x)', [2,5])
            
