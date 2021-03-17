// [1] $ g++ falsa_posicion.cpp -o fp.out -std=c++11 -I/usr/include/python3.8 -lpython3.8 -lginac -lcln
// [2] $ ./fp.out

#include <iostream>
#include <vector>
#include <cmath>
#include <ginac/ginac.h>
#include "matplotlibcpp.h"

using namespace std;
using namespace GiNaC;
namespace plt = matplotlibcpp;

symbol x("x");
ex func;

/**
 * @brief Función no lineal a la que deseamos encontrar solución.
 * @param x_k Variable independiente (Real).
 * @return Retorna el resultado de evaluar f(x).
 */ 
double f(double x_){
    ex x_k = x_;
    ex result = evalf(subs(func,x==x_k));
    return ex_to<numeric>(result).to_double(); 
}


/**
 * @brief Calcula el siguiente valor de la sucesión mediante la fórmula de la secante.
 * @param x_k Valor actual de la sucesión.
 * @param x_prev Valor anterior de la sucesión x_(k-1).
 * @return Retorna el siguiente valor de la sucesión x_(k+1).
 */ 
double calcular_sgte_valor(double x_k, double x_prev){
    return x_k - ((x_k - x_prev)* f(x_k))/(f(x_k) - f(x_prev));
}


/**
 * @brief Calcula el error de la aproximación obtenida. 
 * Entre más cercana a f(x)=0, mejor es la aproximación.
 * @param x_k Valor de la solución apróximada.
 * @return Retorna el resultado en valor absoluto se evaluar f(x_k). 
 */ 
double calcular_error(double x_k){
    return abs(f(x_k));
}


/**
 * @brief Determina si se cumple alguna de las condiciones de parada del método iterativo. 
 *          (1) Se consigue un valor de tolerancia dado.
 *          (2) El denominador se aproxima demasiado a cero.
 * @param x_k Valor actual de la sucesión.
 * @param x_prev Valor anterior de la sucesión x_(k-1).
 * @param tol Valor de tolerancia.
 * @return Verdadero si se cumple alguna condición, falso si no.
 */ 
bool condicion_parada(double x_k, double x_prev, double tol){
    double val = 10e-8;
    if (calcular_error(x_k) < tol) return true;             
    else if (abs(f(x_k) - f(x_prev)) < val) return true;    
    else return false;
}


/**
 * @brief Verifica el teorema de bolzano para la función en un intervalo definido.
 * @param a Valor de inicio del intervalo.
 * @param b Valor final del intervalo.
 * @return Verdadero si se cumple el teorema, falso si no.
 */
bool teorema_bolzano(double a, double b){
    return f(a)*f(b) < 0;
}


/**
 * @brief Grafica el error en funcion de la cantidad de iteraciones.
 * @param x Set de valores del eje x.
 * @param y Set de valores del eje y.
 */ 
void plot(vector<double> x, vector<double> y){
    
    plt::named_plot("Error |f(x_k)|",x,y);
    plt::title("Error Falsa Posición");
    plt::legend();
    plt::show();
}


/**
 * @brief Permite encontrar los ceros de una función f(x) de forma iterativa
 * aplicando el método de la falsa posición.
 * @param a Valor de inicio del intervalo.
 * @param b Valor final del intervalo.
 * @param tol Valor de la tolerancia de resultado aceptable.
 * @param max_itr Cantidad máxima de iteraciones que se pueden realizar.
 */ 
int falsa_posicion(double a, double b, double tol, int max_itr=100){
    
    if (!teorema_bolzano(a,b)){
        throw "[Error 001] La función no cumple con el teorema de Bolzano en el intervalo dado.";
    }

    cout << max_itr << endl;
    vector<double> error, iter;
    double x_k;
    int k = 0;

    while(!condicion_parada(b,a,tol) && k < max_itr){

        x_k = calcular_sgte_valor(b,a);
        error.push_back(calcular_error(x_k));
        iter.push_back(k);

        if (teorema_bolzano(a, x_k)) b = x_k;
        else if (teorema_bolzano(x_k, b)) a = x_k;
        else {
            cout << "[Error 002] No es posible continuar calculando. No se cumple el teorema de Bolzano." << endl;
            return 1;
        }
        k++;
    }

    cout << "Cantidad de iteraciones: " << k << endl;
    cout << "Aproximación de la solución: " << x_k << endl;
    cout << "Error: " << error.at(k-1) << endl;

    plot(iter,error);

    return 0;
}


/**
 * @brief Rutina principal. Obtiene los parámetros de entrada para ejecutar el
 *        método iterativo que aproxima la solución de f(x)=0.
 */ 
void run(){
    symtab table;
    table["x"] = x;
    parser reader(table);

    string ftext;
    cout << "Escriba la función: " << endl;
    cin >> ftext;

    func = reader(ftext);

    double a,b,tol;
    int max_itr;

    cout << "Escriba el valor de a: " << endl;
    cin >> a;

    cout << "Escriba el valor de b: " << endl;
    cin >> b;

    cout << "Tolerancia: " << endl;
    cin >> tol;

    cout << "Cantidad máxima de iteraciones: " << endl;
    cin >> max_itr;

    falsa_posicion(a,b,tol,max_itr);
}


int main(int argc, char const *argv[])
{
    run();
    return 0;
}
