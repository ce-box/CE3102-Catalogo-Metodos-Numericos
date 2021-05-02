/* ------------------------------------------------------------
 * @file: jacobi.cpp
 * @dependencias: armadillo, matplotlibcpp
 * @version 0.1
 * ------------------------------------------------------------*/

// [1] $ g++ jacobi.cpp -o jac.out -std=c++11 -I/usr/include/python3.8 -lpython3.8 
// [2] $ ./jac.out

#include <armadillo>
#include "matplotlibcpp.h"

using namespace std;
using namespace arma;
namespace plt = matplotlibcpp;


/**
 * @brief Indica si una matriz es diagonal dominante
 * @param M Matriz de coeficientes, debe ser cuadrada.
 * @return True si la matriz es diagonal dominante, False si no lo es.
 */
bool dominant_diagonal(mat M){
    int n = M.n_rows;
    double d_value, row_sum;

    for (int i = 0; i < n; i++){
        d_value = M(i,i);
        row_sum = 0; 
        for (int j = 0; j < n; j++){
            row_sum += abs(M(i,j));
        }
        
        if (abs(d_value) < row_sum-abs(d_value))
            return false; 
    }

    return true;
}


/**
 * @brief Retorna la matriz superior
 */
mat upper_mat(mat A){
    int n = A.n_rows;
    mat U(n,n);
    U.zeros();

    for(int i = 0; i < n; i++){
        for (int j = i; j < n; j++){
            if (i == j){
                continue;
            }
            U(i,j) = A(i,j);
        }
    }

    return U;
}


/**
 * @brief Retorna la matriz diagonal
 */
mat diag_mat(mat A){
    int n = A.n_rows;
    mat D(n,n);
    D.zeros();

    for(int i = 0; i < n; i++){
        D(i,i) = A(i,i);
    }

    return D;
}


/**
 * 
 */
mat lower_mat(mat A){
    int n = A.n_rows;
    mat L(n,n);
    L.zeros();

    for(int i = 0; i < n; i++){
        for (int j = 0; j < i; j++){
            L(i,j) = A(i,j);
        }
    }

    return L;
}


/**
 * @brief Calcula la aproximación a la solución de un sistema de ecuaciones por 
 *        el método de Jacobi.
 * @param A Matriz de coeficientes.
 * @param b Vector de términos independientes.
 * @param x_o Vector de valor inicial.
 * @param max_itr Iteraciones máximas.
 */ 
void jacobi(mat A, vec b, vec x_o, int max_itr=100){

    if (!dominant_diagonal(A)){
        cout<<"[Error] La matriz no es diagonal dominante."<<endl;
        return;
    }


}


int main(int argc, char const *argv[])
{
    mat A = {{-10,  3,  0,  0},
             {  4, 15,  2,  1},
             {  3,  8,-21,  4},
             {  1,  1,  1,  4}};

    // Probar Diagonal Dominante
    cout<<"is M DD?: "<< dominant_diagonal(A)<<endl;
    
    // Probar matrices LDU
    mat U = upper_mat(A);
    U.print("U: ");
    mat D = diag_mat(A);
    D.print("D: ");
    mat L = lower_mat(A);
    L.print("L: ");

    return 0;
}
