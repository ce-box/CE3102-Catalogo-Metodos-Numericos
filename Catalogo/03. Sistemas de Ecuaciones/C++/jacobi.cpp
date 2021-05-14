/* ------------------------------------------------------------
 * @file: jacobi.cpp
 * @dependencias: armadillo, matplotlibcpp
 * @version 0.1
 * ------------------------------------------------------------*/

// [1] $ g++ jacobi.cpp -o jac.out -std=c++11 -I/usr/include/python3.8 -lpython3.8 
// [2] $ ./jac.out
// Nota: Asegurarse de tener armadillo bien instalado.

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
 * @brief Retorna la matriz inferior
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
 * @brief Grafica el error en funcion de la cantidad de iteraciones.
 * @param x Set de valores del eje x.
 * @param y Set de valores del eje y.
 */ 
void plot(vector<double> x, vector<double> y){
    
    plt::named_plot("Error |Ax-b|",x,y);
    plt::title("Error Jacobi");
    plt::legend();
    plt::show();
}


/**
 * @brief Calcula la aproximación a la solución de un sistema de ecuaciones por 
 *        el método de Jacobi.
 * @param A Matriz de coeficientes.
 * @param b Vector de términos independientes.
 * @param x_o Vector de valor inicial.
 * @param tol Tolerancia de la aproximación.
 * @param max_itr Iteraciones máximas.
 */ 
void jacobi(mat A, vec b, vec x_o, double tol, int max_itr=100){

    if (!dominant_diagonal(A)){
        cout<<"[Error] La matriz no es diagonal dominante."<<endl;
        return;
    }

    mat L = lower_mat(A);
    mat D = diag_mat(A);
    mat U = upper_mat(A);
    
    vec x = x_o;
    int k = 0;

    mat T = -D.i()*(L+U);
    vec c = D.i()*b;

    vector<double> errors, itr;
    double error = tol;

    while(error > tol && k < max_itr){
        x = T*x + c;
        error = norm(A*x-b);
        errors.push_back(error);
        itr.push_back(k);
        k++;
    }

    x.print("x: ");
    plot(itr,errors);
}


int main(int argc, char const *argv[])
{
    mat M = {{5,1,1},{1,5,1},{1,1,5}};
    vec b = {7,7,7};
    vec x_o = {0,0,0};
    jacobi(M,b,x_o,5);

    return 0;
}
