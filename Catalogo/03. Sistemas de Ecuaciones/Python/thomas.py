import numpy
import numpy as np
def thomas(A, b):
    n = A.shape[0]
    m = A.shape[1]
    t = b.shape[0]
    o = b.shape[1]

    if (n != m):
        raise ValueError("La matriz no es cuadrada")
    if t!=n or o!=1:
        raise ValueError("La matriz b no tiene las dimensiones correctas")

    M=obtiene_verifica_matriz(A,n)
    x=thomas_aux(M[0],M[1],M[2],b,n)

    return x

def obtiene_verifica_matriz(A,n):
    a = np.zeros((n, 1))
    b = np.zeros((n, 1))
    c = np.zeros((n, 1))

    for i in range(0,n):
        for j in range(0,n):
            if j == i:
                a[i] = A[i,j]
            elif (i+1) == j:
                b[i] = A[i,j]
            elif (i-1) == j:
                c[i] = A[i,j]
            else:
                if A[i,j] != 0:
                    raise ValueError("Esta matriz no es tridiagonal")
    return a,b,c

def thomas_aux(a,b,c,d,n):
    r=np.zeros((n, 1))
    t= np.zeros((n, 1))
    sol= np.zeros((n, 1))
    r[0]=b[0]/a[0]

    for i in range(1,n-1):
        if(a[i]-r[i-1]*c[i])==0:
            raise ValueError("No se puede dividir entre 0")
        else:
            r[i]=b[i]/(a[i]-r[i-1]*c[i])
    t[0]=d[0]/a[0]
    for j in range(1,n):
        if (a[j]-r[j-1]*c[j])==0:
            raise ValueError("No se puede dividir entre 0")
        else:
            t[j]=(d[j]-t[j-1]*c[j])/(a[j]-r[j-1]*c[j])

    sol[n-1]=t[n-1]
    k=n-2
    while k>=0:
        sol[k]=t[k]-r[k]*sol[k+1]
        k=k-1
    return sol


if __name__ == '__main__':
    A=np.matrix('5 1 0 0;1 5 1 0;0 1 5 1;0 0 1 5',float)
    b=np.matrix('-12;-14;-14;-12',float)
    print("Solucion x de ejemplo:")
    print(thomas(A, b))