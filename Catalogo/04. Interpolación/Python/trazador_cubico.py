import numpy
import numpy as np
import sys
from sympy import *

#Parametros de entrada
#XK vector de puntos de tamano n, #YK vector de imagenes de tamano n
#Salida
#S vector de polinomios del trazador cubico
#Trazador cubico
def traz_cubico(Xk, Yk):
    n = len(Xk) #Largo del vector

    if (len(Xk) != len(Yk)): #Comprueba que los vectores sean iguales
        return ('Los venctores no cumplen con la condicion de tamano')
    h = np.zeros(n-1)

    for i in range(0, n-1): #Calculo del vector h (distancia entre cada punto
                            # del trazador)
        h[i] = Xk[i + 1] - Xk[i]

    A = np.zeros((n-2,n-2)) #Matriz tridiagonal

    A[0][0] = 2*(h[0]+h[1]) #Calculo de la primer fila de la matriz n-1xn-1
    A[0][1] = h[1]

    for i in range(1, n-3): #Calculo del la de 1 a n-2 de la matriz
        A[i][i] = h[i]
        A[i][i+1] = 2*(h[i] + h[i+1])
        A[i][i+2] = h[i+1]
    
    A[n-3][n-4] = h[n-3]  #Calculo de la posicion n-1 de la matriz
    A[n-3][n-3] = 2*(h[n-3] + h[n-2])

    u = np.zeros(n-2) # Vector de variables u

    for i in range(0, n-2): #Calculo del vector u
        u[i] = 6*(((Yk[i+2]-Yk[i+1])/h[i+1])-((Yk[i+1]-Yk[i])/h[i]))

    M_temp = thomas(A,u) #Implementacion del metodo de Thomas para resolver
                         #El sistema Ax=U
    M = np.zeros(n)
    M[0] = 0 #Matriz M con M0 = 0 y Mn-1 = 0
    M[n-1] = 0

    for i in range(1, n-1): #Construccion de la matriz M apartir de la matriz 
                            #M_temp
        M[i] = M_temp[i-1]
    
    S = [] #Vector de polinimios de solucion S
    
    for i in range(0, n-1): #Calculo de las variables a,b,c,d
        a = (M[i+1]-M[i])/(6*h[i])
        b = M[i]/2
        c = (Yk[i+1]-Yk[i])/h[i] - (h[i]/6)*(M[i+1]+2*M[i])
        d = Yk[i]
        S.append(crear_funcion(a,b,c,d,Xk[i])) #Creacion de polinomio Si 

    print(S)
    return(S)

#Parametros de entrada
# a,b,c,d valores del polinimio
# xk valor de x0
#Salida
#polinomio del trazador cubico Si
#Estructura del polinomio Si   
def crear_funcion(a,b,c,d,xk):
    x = Symbol('x')
    polinomio = 0
    polinomio = sympify(polinomio)
    polinomio = a*(x - xk)**3 + b*(x - xk)**2 + c*(x - xk) + d 
    return (polinomio)

#Parametros de entrada
#A matriz nxn, #b matriz de valores independientes
#Salida
#x matriz  de las soluciones
#Metodo de thomas
def thomas(A, b):
    n = A.shape[0]  # Obtiene el numero de filas de A
    m = A.shape[1]  # Obtiene el numero de columnas de A
    
    M=obtiene_verifica_matriz(A,n)#Metodo para comprobar si la matriz es tridiagonal o lanza error si no lo es
    x=thomas_aux(M[0],M[1],M[2],b,n)#Llama a la funcion auxiliar para resolver el problema
    return x
#Parametros de entrada
#A matriz nxn, #n largo de la matriz
#Salida
#a Matriz de los valores de la diagonal,#b Matriz con los valores encima de la diagonal,#c Matriz con los valores debajo de la diagonal
#Metodo de thomas
def obtiene_verifica_matriz(A,n):
    a = np.zeros((n, 1))#Matriz lleno de zeros nx1
    b = np.zeros((n, 1))
    c = np.zeros((n, 1))

    for i in range(0,n):#Mueve filas
        for j in range(0,n):#Mueve columnas
            if j == i:#Diagonal
                a[i] = A[i,j]#Guarda valor en la matriz
            elif (i+1) == j:
                b[i] = A[i,j]
            elif (i-1) == j:
                c[i] = A[i,j]
            else:
                if A[i,j] != 0:
                    raise ValueError("Esta matriz no es tridiagonal")
    return a,b,c


#Parametros de entrada
#a Matriz de los valores de la diagonal,#b Matriz con los valores encima de la diagonal,#c Matriz con los valores debajo de la diagonal
#Salida
#sol matriz de soluciones
#Metodo de auxiliar de thomas, realiza el algoritmo para resolver el sistema
def thomas_aux(a,b,c,d,n):
    r=np.zeros((n, 1))#Matriz lleno de zeros nx1
    t= np.zeros((n, 1))
    sol= np.zeros((n, 1))
    r[0]=b[0]/a[0]#Primer coeficiente

    for i in range(1,n-1):
        if(a[i]-r[i-1]*c[i])==0:#Comprueba que el divisor no sea 0
            raise ValueError("No se puede dividir entre 0")
        else:
            r[i]=b[i]/(a[i]-r[i-1]*c[i])#Calcula los nuevos coeficientes
    t[0]=d[0]/a[0]#Se realiza el primera barrido/susticion(similar hacia adelante)
    for j in range(1,n):
        if (a[j]-r[j-1]*c[j])==0:#Comprueba que el divisor no sea 0
            raise ValueError("No se puede dividir entre 0")
        else:
            t[j]=(d[j]-t[j-1]*c[j])/(a[j]-r[j-1]*c[j])#Completa el barrido/susticion(similar hacia adelante)

    sol[n-1]=t[n-1]#Calcula la ultima solucion
    k=n-2
    while k>=0:
        sol[k]=t[k]-r[k]*sol[k+1]#Calcula las demas soluciones
        k=k-1
    return sol

# Valores iniciales
if __name__ == '__main__':
    Xk = [1,1.05,1.07,1.1]
    Yk = [2.718282, 3.286299, 3.527609, 3.905416]
    traz_cubico(Xk, Yk)