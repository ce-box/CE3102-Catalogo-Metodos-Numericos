/* ------------------------------------------------------------
 * @file: fact_cholesky.cpp
 * @dependencias: None
 * @version 3.0
 * ------------------------------------------------------------*/

// Compile and run:
// [1] $ g++ fact_cholesky.cpp -o fc.out -std=c++11 
// [2] $ ./fc.out

#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;


/**
 * @brief Determina si una matriz es cuadrada, es decir, si posee
 *        el mismo número de columnas y filas.
 * @param A Matriz a ser verificada.
 * @return True si la matriz es cuadrada, false si no lo es.
 */ 
bool is_square(mat A){
    return A.n_rows == A.n_cols;
}


/**
 * @brief Verifica que una matriz sea simétrica.
 * @param A Matriz de tamaño nxn.
 * @return True si es simétrica, false si no lo es.
 */ 
bool is_simetric(mat A){
    return approx_equal(A, A.t(),"absdiff", 0);
}


// Source: https://www.tutorialspoint.com/cplusplus-program-to-compute-determinant-of-a-matrix
// Nota: Se implementó porque la función de armadillo arma::det(mat) presentó problemas.
/**
 * @brief Calcula el determinante de una matriz cuadrada.
 * @param A Matriz de tamaño nxn.
 * @param n Tamaño de la matriz.
 * @return El valor del determinante de la matriz.
 */ 
double determinant(mat A, int n) {
    double det = 0;
    mat sub_matrix(10,10);
    
    if (n == 1){
        // Caso 0: Si la matriz es de orden 1
        return A(0,0);
    } else if (n == 2){
        // Caso base: Determinante de una matriz 2x2
        return (A(0,0) * A(1,1) - A(1,0) * A(0,1));
    } else {
        // Caso de reducciones de la matriz
        for (int x = 0; x < n; x++){
            int sub_i = 0;
            for (int i = 1; i < n; i++){
                int sub_j = 0;
                for (int j = 0; j < n; j++){
                    if (j == x)
                        continue;
                    sub_matrix(sub_i,sub_j) = A(i,j);
                    sub_j++;
                }
                sub_i++;
            }
            det += (pow(-1,x)* A(0,x) * determinant(sub_matrix, n-1));
        }
    }
    return det;
}


/**
 * @brief Determina si una matrix dada es Definida Positiva.
 * @param A Matriz de tamaño nxn.
 * @return True si es definida positiva, false si no es.
 */ 
bool is_positive_definite(mat A){
    int n = A.n_rows;

    for (int i = 0; i < n; i++){
        mat A_t = A.submat(0,0,i,i);
        double det_result = determinant(A_t,i+1);
        if (det_result < 0)
            return false;
    }
    return true;
}



/**
 * @brief Resuelve un sistema de ecuaciones por el método de la sustitución
 *        hacia atrás.
 * @param A Matriz de tamaño nxn. Debe ser triangular superior.
 * @param b Vector de tamaño n.
 * @return El vector con la solución al sistema de ecuaciones.
 */
vec back_substitution(mat A, vec b){
    int n = A.n_rows;
    vec x(n);
    x.zeros();

    for (int i = n-1; i >= 0; i--){
        double sum = 0;
        for (int j = i+1; j < n; j++){
            sum += A(i,j) * x(j);
        }
        x(i) = (1/A(i,i)) * (b(i) - sum);
    } 

    return x;
}


/**
 * @brief Resuelve un sistema de ecuaciones por el método de la sustitución
 *        hacia adelante.
 * @param A Matriz de tamaño nxn. Debe ser triangular inferior.
 * @param b Vector de tamaño n.
 * @return El vector con la solución al sistema de ecuaciones.
 */
vec forward_substitution(mat A, vec b){
    int n = A.n_rows;
    vec x(n);
    x.zeros();

    for (int i = 0; i < n; i++){
        double sum = 0;
        for (int j = 0; j < i; j++){
            sum += A(i,j) * x(j);
        }
        x(i) = (1/A(i,i)) * (b(i) - sum);
    } 

    return x;
}


/**
 * @brief Soluciona un sistema de ecuaciones utilizando la Factorización de Cholesky.
 * @param A Matriz de coeficientes. Deber ser cuadrada, simétrica y positiva definida.
 * @param b Vector de términos independientes.
 */ 
void fact_cholesky(mat A, vec b){
    
}




int main(int argc, char const *argv[]){

    // Test sustitución hacia atrás
    mat A(4,4);
    A = {{1, 1,-1, 3},
         {0,-1,-1,-5},
         {0, 0, 3,13},
         {0, 0, 0,-13}};
    vec b(4);
    b = {4, -7, 13, -13};
    vec x = back_substitution(A,b);
    x.print("x(bs): ");


    // Test sustitución hacia adelante
    mat C(5,5);
    vec d(5);
    C = {{2, 0, 0, 0, 0},
         {3, 2, 0, 0, 0},
         {3,-1,-1, 0, 0},
         {1, 1, 1, 1, 0},
         {1, 2, 3,-4, 5}};
    d = {-8, 10 , 2, 1, 3};
    x = forward_substitution(C,d);
    x.print("x(fs): ");

    // Test simetrica
    mat E(3,3);
    E = {{4,2,1},
         {2,5,2},
         {1,2,6}};
    cout<<"is E simetric?: "<< is_simetric(E)<< endl;

    // Test diagonal positiva
    cout<<"Det(E): "<< determinant(E,3)<< endl;
    cout<<"is E DP?: "<< is_positive_definite(E)<< endl;

    return 0;
}
