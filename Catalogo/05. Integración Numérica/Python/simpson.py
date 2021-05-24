import numpy as np
import sys
from sympy import *

#Parametros de entrada
#func: funcion integrable, v: interbalo del funcion 
#Salida
#f_resultado: aproximación de la integral, error: error de la aproximacion
#Funcion de Simpson
def simpson(func, v):

    a = v[0] #Constantes
    b = v[1]

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy  
    fx =  lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)

    x0 = a
    x1 = (a+b)/2
    x2 = b

    h = (b-a)/2 #Valor h

    f_resultado = (h/3)*(fx(x0) + 4*fx(x1) + fx(x2)) #Resultado de la aproximacion
    
    df = f.diff(x,4) #Cuarta derivada 
    dfx = lambdify(x, df, modules=['numpy']) #Se iniacilaiza la funcion f4(x)

    max_relativo = maximo(dfx, v) #Calculo del maximo en los extremos del intervalo

    #En casa de que la funcion no se indefina f(x) = 0
    try:
        sol = np.solve(f.diff(x,1))
        punto = []
        for i in sol: #Calcula cual de los resultados es el maximo
            if (abs(dfx(i)) > punto):
                punto = [i, abs(dfx(i))]
        
        if (punto[1] > max_relativo[1]): #Compara el maximo con el maximo de los extremos
            error = (h**5/90)*abs(dfx(punto[0]))
        else:
            error = (h**5/90)*abs(dfx(max_relativo[0])) 
        
    except:
        error = (h**5/90)*abs(dfx(max_relativo[0]))
    
    print(f_resultado, error)
    return(f_resultado, error) #Resultado (aproximacion, error)

#Parametros de entrada
#f: cuarta derivada de la funcion, v: interbalo del funcion 
#Salida
#[v,f*] punto maximo, valor del punto maximo 
#Calculo del maximo entre dos puntos de una funcion
def maximo(f, v):
    f1 = abs(f(v[0])) #Valores de los puntos
    f2 = abs(f(v[1]))

    if (f1 > f2): #Punto maximo 
        return [v[0], f1] #Retorno [punto, maximo]
    else:
        return [v[1], f2] #Retorno [punto, maximo]

#Valores iniciales 
if __name__ == '__main__':
    f = 'ln(x)' #Funcion
    v = [2,5]   #Intervalo 
    simpson(f, v)