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
    return 0;
}
