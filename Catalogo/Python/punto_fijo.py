import numpy as np
import sys
from sympy import *
from matplotlib import pyplot as plt

def punto_fijo():

    #Funcion de prueba: "E**x - x - 2", [0,2], 10**(-10), 100

    x = Symbol('x') #Inicializa "x" como símbolo
    f = sympify(func) #Se traduce el string "func" a una función de sympy 
    fx = lambdify(x, f, modules=['numpy']) #Se inicializa la función f(x)

    D = []
    A = []

    cont = 0 #Se inicializa el contador

    
