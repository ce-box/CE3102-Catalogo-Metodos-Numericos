/* ------------------------------------------------------------
 * @file: pseudoinversa.cpp
 * @dependencias: armadillo, matplotlibcpp
 * @version 0.1
 * ------------------------------------------------------------*/

// [1] $ g++ pseudoinversa.cpp -o psi.out -std=c++11
// [2] $ ./psi.out
// Nota: Asegurarse de tener armadillo bien instalado.

#include <armadillo>

using namespace std;
using namespace arma;


/**
 * @brief Calcula la aproximación a la solución de un sistema de ecuaciones por 
 *        el método de la pseudoinversa.
 * @param A Matriz de coeficientes.
 * @param b Vector de términos independientes.
 * @param tol Tolerancia de la aproximación.
 * @param max_itr Iteraciones máximas.
 */
void pseudoinversa(mat A, vec b, double tol, int max_itr=15){
    int n = A.n_rows;
    int m = A.n_cols;
    double alpha = eig_sym(A*A.t()).max();

    mat I(n,m);
    I.eye();
    mat x = (1/alpha)*A.t();

    double error = tol;
    int k = 0;

    while (k < max_itr){
        x = x*(2*I-A*x);
        error = norm((A*x*A)-A);

        if (error < tol)
            break;
        k++;
    }
    mat x_p = x*b;

    // Mostrar los resultados
    x_p.print("x_p: \n");
    cout<<"Error: "<< norm(A*x_p-b)<<endl;

}


int main(int argc, char const *argv[])
{
    mat A = {{ 1, 2, 4},
             { 2,-1, 1},
             { 1, 0, 1}};
    vec b = {4,3,9};
    pseudoinversa(A,b,10e-8,3);
    return 0;
}
