#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>

using namespace std;

/**
 * @brief Función no lineal a la que deseamos encontrar solución.
 * @param x Variable independiente (Real).
 * @return Retorna el resultado de evaluar f(x).
 */ 
double f(double x){
    // f(x) = cos²(2x) - x²
    return pow(cos(2*x), 2) - pow(x,2); 
}


/**
 * @brief Calcula el siguiente valor de la sucesión.
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
 * 
 */ 
bool condicion_parada(double x_k, double x_prev, double tol){
    if (calcular_error(x_k) < tol) return true;
    else if ((f(x_k) - f(x_prev)) < 10e-8) return true;
    else return false;
}


/**
 * 
 */ 
tuple<double, vector<double>> secante(double x_0, double x_1, double tol, int max_itr=100){

    vector<double> error;
    double x_k = x_1, x_prev = x_0, x_sgt;
    int k = 1;

    while (!condicion_parada(x_k, x_prev, tol) || k > max_itr){
        cout << "x_k = " << x_k << endl;
        cout << "x_k-1 = " << x_prev << endl;
        x_sgt = calcular_sgte_valor(x_k, x_prev);
        x_prev = x_k;
        x_k = x_sgt;
        error.push_back(calcular_error(x_k));
        k++;
    }

    return make_tuple(x_k, error);
}



int main(int argc, char const *argv[])
{
    tuple<double, vector<double>> r = secante(0.75, 0.5, 0.00001);

    // Resultado
    cout << "aprox:: " << get<0>(r) << endl;

    return 0;
}
