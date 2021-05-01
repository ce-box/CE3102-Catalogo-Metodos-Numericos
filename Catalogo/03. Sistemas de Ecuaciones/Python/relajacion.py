import numpy as np
from sympy import Symbol
import math


def relajacion(matriz, vector_independiente, valor_inicial, tol, iterMax):
    """
    Funcion que implementa el metodo iterativo de relajacion
    :param matriz_a: Matriz cuadrada sim´etrica definida positiva
    :param vector_independiente: Vector columna
    :param tol: Tolerancia al fallo que debe tener el resultado
    :param iterMax: Iteraciones maximas
    :return: Vector resultado
    """
    n = len(matriz)  # Numero de filas
    m = len(matriz[0])  # Numero de columnas
    iter = 0  # Numero de iteraciones
    error = 0  # Error de la aproximacion

    # Verificar si es cuadrada
    if (n != m):
        return 'La matriz debe ser cuadrada'

    # Calcular la matriz D
    d = cal_matriz_d(matriz, n)

    # Calcular la matriz L
    l = cal_matriz_l(matriz, n)

    # Calcular matriz U
    u = cal_matriz_u(matriz, n)

    # Calcular N
    matriz_n = -(l + u)

    # Calcular M-1
    m_inver = np.linalg.inv(d)

    # Calcular valor M-1*N
    matriz_a = np.dot(m_inver, matriz_n)

    # Calcular radio espectral
    p = radio_espectral(matriz_a)

    # Calcular el valor de w optimo
    w = 2 / (1 + math.sqrt(1 - p ** 2))

    # Calcula d+wl
    temp_1 = d + np.dot(w, l)

    # Calcula w*b
    temp_3 = np.dot(w, vector_independiente)

    # Calculo de C
    C = hacia_Adelante(temp_1, temp_3, n)


    while (iter <= iterMax):

        if (error > tol):
            break

        temp_2 = np.dot((np.dot((1 - w), d) - np.dot(w, u)), valor_inicial)

        T = hacia_Adelante(temp_1, temp_2, n)  # Se calcula T
        X = T + C
        iter += 1
        valor_inicial = X
        print(valor_inicial)


# Parametros de entrada
# matriz_a matriz  cuadrada, #n largo de la matriz
# Salida
# Matriz D
# Calculo de la matriz D
def cal_matriz_d(matriz_a, n):
    d = np.zeros_like(matriz_a)
    for i in range(0, n):
        d[i][i] = matriz_a[i][i]
    return d


# Parametros de entrada
# matriz_a matriz  cuadrada, #n largo de la matriz
# Salida
# Matriz L
# Calculo de la matriz L
def cal_matriz_l(matriz_a, n):
    d = np.zeros_like(matriz_a)
    for i in range(0, n):
        for j in range(0, i):
            d[i][j] = matriz_a[i][j]
    return d


# Parametros de entrada
# matriz_a matriz  cuadrada, #n largo de la matriz
# Salida
# Matriz U
# Calculo de la matriz U
def cal_matriz_u(matriz_a, n):
    u = np.zeros_like(matriz_a)
    for i in range(0, n):
        for j in range(i + 1, n):
            u[i][j] = matriz_a[i][j]
    return u


# Parametros de entrada
# matriz_a matriz  cuadrada, matriz cuadrada D
# Salida
# Matriz N
# Calculo de la matriz N
def cal_matriz_n(matriz_a, matriz_d):
    matriz_n = -(matriz_a - matriz_d)
    return matriz_n


# Parametros de entrada
# matriz_a matriz  cuadrada
# Salida
# Valor del radio espectral
# Calculo del valor máximo del radio espectral
def radio_espectral(matriz_a):
    valores, vectores = np.linalg.eig(matriz_a)
    mayor = 0
    for i in range(len(valores)):
        if (mayor < abs(valores[i])):
            mayor = valores[i]
    return mayor


#Parametros de entrada
#A matriz  triangular inferior, #b nueva matriz de valores independientes, #n largo de la matriz
#Salida
#soluciones matriz de soluciones "x#"
#Sustitucion hacia adelante
def hacia_Adelante(A, b, n):
    soluciones = np.zeros((n, 1))#Matriz de soluciones
    soluciones[0] = b[0] / A[0, 0]#Primera solucion
    i = 1

    while (i < n):#Mueve las filas
        sumatoria = 0
        for j in range(0, i):#Mueve columnas
            sumatoria = sumatoria + (A[i, j] * soluciones[j])#Sumatoria necesaria

        x = (b[i] - sumatoria) / A[i, i]#Calcular el valor de la solucion "x"
        soluciones[i] = x#Agrega la solucion a la matriz de soluciones
        i = i + 1
    return soluciones

if __name__ == '__main__':
    a = [[4, 3, 0], [3, 4, -1], [0, -1, 4]]
    b = [[7], [7], [-1]]
    relajacion(a, b, [[0], [0], [0]], 10, 5)