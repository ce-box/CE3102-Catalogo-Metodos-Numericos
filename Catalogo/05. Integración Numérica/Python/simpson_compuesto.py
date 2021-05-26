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

    a = interv[0] #Se extrae el valor inicial del intervalo
    b = interv[-1] #Se extrae el valor final del intervalo
    m = 7 #Cantidad de puntos

    h = (b-a)/(m-1) #Se calcula el valor de "h"
    x = [] #Se inicializa la lista de los valores de x (x0, x1, x2, ...)
    
    cont = 0 #Inicializa contador para bucle

    while(cont < m):
        
        xk = a+cont*h #Se calcula xk
        x.append(xk) #Se añade xk a la lista de valores de x
        cont += 1 #Aumenta contador

    cont = 1 #Inicializa contador para bucle (en 1 porque no contempla x0)
    sumPar = 0 #Inicializa resultado de sumatoria de indices pares
    sumImp = 0 #Inicializa resultado de sumatoria de indices impares
            
    while(cont < (len(x) - 1)): #Punto de parada "len(x) - 1" porque no contempla xn

        if(par(cont)): #Verifica indice par
            sumPar += fx(x[cont]) #Sumatoria indices pares
        else: #Si el indice no es par
            sumImp += fx(x[cont]) #Sumatoria indices impares
        cont += 1 #Incrementa contador

    aprox = (h/3)*(fx(x[0]) + 2*sumPar + 4*sumImp + fx(x[-1])) #Calcula aproximacion
    error = (((b - a)*(h**4))/180) * abs(d4fx(2)) #Calcula error

    print([aprox, error]) #Muestra en pantalla la aproximacion y el error
    return[aprox, error] #Retorna aproximacion y error
    

def par(indice): #Verifica paridad de indice
    if((indice%2) == 0): #Si la division entre 2 da cero (es par)
        return True #True si es par
    else: #Si no da cero (es impar)
        return False #False si es impar

simpson_compuesto('log(x)', [2,5]) #Ejemplo
            
