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



void trapecio_compuesto(ex f, double a, double b, int m){

    symbol x("x");
    double h = (b - a)/(m-1);
    
    double x_prev = a, x_n;
    double sum = 0;

    for (int n = 0; n < m; n++ ){
        x_n = x_prev + n*h;
        sum += f.subs(x == x_prev) + f.subs(x == x_n);
        
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
    /* code */
    return 0;
}
