/* ------------------------------------------------------------
 * @file: trapecio_compuesto.cpp
 * @dependencias: GiNaC
 * @version 0.1
 * ------------------------------------------------------------*/

// [1] $ g++ trapecio_compuesto.cpp -o trp_cmp.out -std=c++11 -I/usr/include/python3.8 -lpython3.8 -lginac -lcln
// [2] $ ./trp_cmp.out

#include <iostream>
#include <vector>
#include <cmath>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

symbol x("x");
ex func;

/**
 * @brief Evalua la función en un valor x
 * 
 * @param x_ Variable independiente (Real).
 * @return Retorna el resultado de evaluar f(x).
 */ 
double f(double x_){
    ex x_k = x_;
    ex result = evalf(subs(func,x==x_k));
    return ex_to<numeric>(result).to_double(); 
}


/**
 * @brief 
 * 
 * @param fx 
 * @param a 
 * @param b 
 * @param m 
 */
void trapecio_compuesto(ex fx, double a, double b, int m){

    func = fx;
    double h = (b - a)/(m-1);
    
    double x_o = a; 
    double x_prev = x_o;
    double x_n;
    
    double sum = 0;

    for (int n = 1; n < m; n++ ){
        x_n = x_o + n*h;
        sum += f(x_prev) + f(x_n);
        x_prev = x_n;
    }

    double result = h/2 * sum;

    cout << " I = " << result << endl;
}


double cota_error(){
    return 0;
}


int main(int argc, char const *argv[])
{
    ex f = log(x);
    trapecio_compuesto(f,2,5,4);
    return 0;
}
