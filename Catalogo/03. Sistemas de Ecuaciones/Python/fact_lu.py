import numpy
import numpy as np


def fact_lu(A, b):
    n = A.shape[0]
    m = A.shape[1]
    t = b.shape[0]
    o = b.shape[1]
    if (n != m):
        raise ValueError("La matriz no es cuadrada")
    if t!=n or o!=1:
        raise ValueError("La matriz b no tiene las dimensiones correctas")
    if(teorema_2_lu(A,n)==False):
        raise ValueError("La matriz no cumple el teorema 2 de Fact LU")
    M = matriz_inf_lu(A, n)

    y = hacia_Adelante(M[1], b, n)
    x = hacia_Atras(M[0], y, n)
    return x

def teorema_2_lu(A, n):
    Valor = True
    for k in range(0, n):
        if np.linalg.det(A[0:k, 0:k]) == 0:
            Valor = False
            break
    return Valor

def matriz_inf_lu(A,n):
    U=A
    L=np.zeros((n,n))

    for k in range(0,n-1):
        for i in range(k+1,n):
            Mik=U[i,k]/U[k,k]
            L[i,k]=Mik
            for j in range(k,n):
                U[i,j]=U[i,j]-Mik*U[k,j]
                if i==j:
                    L[i,j]=1
    L[0,0]=1
    return U,L



def hacia_Atras(A, b, n):
    soluciones = np.zeros((n, 1))
    i = n - 1;

    while (i >= 0):
        sumatoria = 0
        for j in range(i + 1, n):
            sumatoria = sumatoria + (A[i, j] * soluciones[j])

        x = (1 / A[i, i]) * (b[i] - sumatoria)
        soluciones[i] = x
        i = i - 1
    return soluciones


def hacia_Adelante(A, b, n):
    soluciones = np.zeros((n, 1))
    soluciones[0] = b[0] / A[0, 0]
    i = 1;

    while (i < n):
        sumatoria = 0
        for j in range(0, i):
            sumatoria = sumatoria + (A[i, j] * soluciones[j])

        x = (b[i] - sumatoria) / A[i, i]
        soluciones[i] = x
        i = i + 1
    return soluciones



if __name__ == '__main__':
    A = np.matrix('4 -2 1; 20 -7 12;-8 13 17', float)
    A2 = np.matrix('25 15 -5 -10;15 10 1 -7;-5 1 21 4;-10 -7 4 18',float)

    b = numpy.matrix('11;70;17')
    b2 = np.matrix('-25;-19;-21;-5')
    print("Ejemplo 1 tiene como resultado:")
    print(fact_lu(A, b))
    print("Ejemplo 2 tiene como resultado:")
    print(fact_lu(A2, b2))

