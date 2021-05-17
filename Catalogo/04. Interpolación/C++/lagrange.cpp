/* ------------------------------------------------------------
 * @file: lagrange.cpp
 * @dependencias: GiNaC
 * @version 0.1
 * ------------------------------------------------------------*/

// [1] $ g++ lagrange.cpp -o lg.out -std=c++11 -I/usr/include/python3.8 -lpython3.8 -lginac -lcln
// [2] $ ./lg.out

#include <iostream>
#include <vector>
#include <cmath>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;


/**
 * @brief Calcula el polinomio Lk (Núcleo de Lagrange).
 * 
 * @param xv Vector de valores de x.
 * @param k Posición del arreglo que se debe saltar para evitar que el denominador sea cero.
 * @return ex Polinomio Lk.
 */
ex lagrange_basis(vector<double> xv, int k){
    
    symbol x("x");
    int n = xv.size();
    ex Lk = 1;

    for (int i = 0; i < n; i++){
        if (i != k){
            Lk *= (x - xv[i])/(xv[k] - xv[i]);
        }
    }

    return Lk;
}


/**
 * @brief Calcula el polinomio de interpolación por el método de Lagrange.
 * 
 * @param xv Vector con valores de x.
 * @param yv Vector con valores de y.
 * @return ex Polinomio de interpolación.
 */
ex lagrange(vector<double> xv, vector<double> yv){

    symbol x("x");
    int n = xv.size();
    ex poly = 0;

    for (int k = 0; k < n; k++){
        poly += yv[k] * lagrange_basis(xv,k);
    }

    return poly.normal();
}


int main(int argc, char const *argv[])
{
    vector<double> xv{-2, 0, 1};
    vector<double> yv{0, 1, -1};
    symbol x("x");
    ex p = lagrange(xv,yv);
    cout << "P(x) = " << p << endl; 
    return 0;
}
