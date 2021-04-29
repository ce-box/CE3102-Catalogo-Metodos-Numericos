/* ------------------------------------------------------------
 * @file: fact_cholesky.cpp
 * @dependencias: None
 * @version 3.0
 * ------------------------------------------------------------*/

// Compile and run:
// [1] $ g++ fact_cholesky.cpp -o fc.out -std=c++11 
// [2] $ ./fc.out

#include <armadillo>

using namespace std;
using namespace arma;


void fact_cholesky(mat A, vec b){
    
}


bool is_simetric(mat A){
    return true;
}


bool is_positive_dependant(mat A){
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
    x.fill(0);

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
    x.fill(0);

    for (int i = 0; i < n; i++){
        double sum = 0;
        for (int j = 0; j < i; j++){
            sum += A(i,j) * x(j);
        }
        x(i) = (1/A(i,i)) * (b(i) - sum);
    } 

    return x;
}


int main(int argc, char const *argv[])
{
    /* code */
    return 0;
}
