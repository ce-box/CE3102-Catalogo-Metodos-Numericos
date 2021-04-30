/* ------------------------------------------------------------
 * @file: fact_cholesky.cpp
 * @dependencias: armadillo, math
 * @version 3.6
 * ------------------------------------------------------------*/

// Compile and run:
// [1] $ g++ fact_cholesky.cpp -o fc.out -std=c++11 
// [2] $ ./fc.out

#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;


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
        //cout << "Det(A_t): "<< det_result<< endl;
        if (det_result < 0){
            return false;
        }
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

    if (is_simetric(A) && is_positive_definite(A)){
        int n = A.n_rows;
        mat L(n,n);
        L.zeros();

        // Calcular la matriz L
        for (int i = 0; i < n; i++){

            for (int j = 0; j < i+1; j++){
                double sum = 0;

                if (i == j){
                    for (int k = 0; k < j; k++)
                        sum += pow(L(j,k),2);
                    L(j,j) = sqrt(A(j,j) - sum);
                } else {
                    for (int k = 0; k < j; k++)
                        sum += L(i,k) * L(j,k);
                    L(i,j) = (A(i,j) - sum) / L(j,j);
                }
            }
        }

        // Resolve Ly=b, para luego solucionar L_t x = y.
        vec y = forward_substitution(L,b);
        vec x = back_substitution(L.t(),y);

        // Presentar la solución del Sistema.
        x.print("x: ");

    } else
        cout<< "[Error 504] El sistema no cumple  con las condiciones para resolverse por Cholesky."<< endl;
}


int main(int argc, char const *argv[]){

    // Test Cholesky
    // ----------------------------
    mat A(4,4);
    vec b(4);

    A = {{ 25, 15, -5,-10},
         { 15, 10,  1, -7},
         { -5,  1, 21,  4},
         {-10, -7,  4, 18}};

    b = {-25, -19, -21, -5};

    cout << "Factorización Cholesky" << endl;
    A.print("A:\n");
    b.print("b:\n");
    fact_cholesky(A,b);

    // Solución del ejemplo -> x = [-1, -1, -1, -1]

    return 0;
}
