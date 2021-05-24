import numpy as np
import sys
from sympy import *

def simpson(func, v):

    a = v[0]
    b = v[1]

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy  
    fx =  lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)

    x0 = a
    x1 = (a+b)/2
    x2 = b

    h = (b-a)/2

    f_resultado = (h/3)*(fx(x0) + 4*fx(x1) + fx(x2))
    
    df = f.diff(x,4)
    dfx = lambdify(x, df, modules=['numpy'])

    max_relativo = maximo(dfx, v)

    try:
        sol = np.solve(f.diff(x,1))
        punto = []
        for i in sol:
            if (abs(dfx(i)) > punto):
                punto = [i, abs(dfx(i))]
        
        if (punto[1] > max_relativo[1]):
            error = (h**5/90)*abs(dfx(punto[0]))
        else:
            error = (h**5/90)*abs(dfx(max_relativo[0]))
        
    except:
        error = (h**5/90)*abs(dfx(max_relativo[0]))
    
    print(f_resultado, error)
    return(f_resultado, error)

def maximo(f, v):
    f1 = abs(f(v[0]))
    f2 = abs(f(v[1]))

    if (f1 > f2):
        return [v[0], f1]
    else:
        return [v[1], f2]

if __name__ == '__main__':
    f = 'ln(x)'
    v = [2,5]
    simpson(f, v)