/* ------------------------------------------------------------
 * @file: trapecio_compuesto.cpp
 * @dependencias: GiNaC
 * @version 2.1
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
 * @brief Evalua la función en un valor x
 * 
 * @param fx Función que se evaluará.
 * @param x_ Variable independiente (Real).
 * @return Retorna el resultado de evaluar f(x).
 */ 
double f(double x_, ex fx){
    ex x_k = x_;
    ex result = evalf(subs(fx,x==x_k));
    return ex_to<numeric>(result).to_double(); 
}


/**
 * @brief Calcula la cota de error de la regla del trapecio compuesto
 * 
 * @param a  Límite inferior.
 * @param b  Límite superior.
 * @param h  Intervalo entre puntos.
 * @return El valor de la cota de error. 
 */
double cota_error(double a, double b, double h){
    
    ex d2_fx = abs(func.diff(x,2));
    double d2_fx_num;

    // Calcular el punto máximo de la función
    double d_fa = f(a,d2_fx);
    double d_fb = f(b, d2_fx);

    if (d_fa > d_fb){
        d2_fx_num = d_fa;
    } else {
        d2_fx_num = d_fb;
    }

    double error = ((b-a) * std::pow(h,2) / 12)* d2_fx_num;
    return error;
}


/**
 * @brief Regla compuesta del trapecio para la calcular la integra de una función.
 * 
 * @param fx Función a integrar con variable x.
 * @param a  Límite inferior.
 * @param b  Límite superior.
 * @param m  Cantidad de puntos a utilizar
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
    double error = cota_error(a,b,h);

    cout << "f(x) = " << func << endl;
    cout << "I = " << result << endl;
    cout << "Error: " << error << endl;

}


// Ejemplo. P10 [Ejemplo 36/47]
int main(int argc, char const *argv[])
{
    ex f = log(x);
    trapecio_compuesto(f,2,5,4);
    
    return 0;
}
