#include <iostream>
#include <vector>
#include <tuple>

using namespace std;

double f(double x){ return 0; }


/**
 * 
 */ 
tuple<double, vector<double>> secante(double x_0, double x_1, double tol, int max_itr=100){

    vector<double> error;
    double x_k = x_1, x_prev = x_0, x_sgt;
    int k = 1;

    while (!condicion_parada(x_k, x_prev, tol) || k < max_itr){
        x_sgt = calcular_sgte_valor(x_k, x_prev);
        x_prev = x_k;
        x_k = x_sgt;
        error.push_back(calcular_error(x_k));
        k++;
    }

    return make_tuple(x_k, error);
}


/**
 * 
 */ 
bool condicion_parada(double x_k, double x_prev, double tol){
    if (calcular_error(x_k) < tol) return true;
    else if (f(x_k) - f(x_prev) < 10e-8) return true;
    else return false;
}


/**
 * 
 */ 
double calcular_sgte_valor(double x_k, double x_prev){
    return x_k - ((x_k - x_prev)* f(x_k))/(f(x_k) - f(x_prev));
}


/**
 * 
 */ 
double calcular_error(double x_k){
    return abs(f(x_k));
}






int main(int argc, char const *argv[])
{
    /* code */
    return 0;
}
