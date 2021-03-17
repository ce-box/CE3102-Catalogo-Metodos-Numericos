#include <iostream>
#include <ginac/ginac.h>

using namespace std;
using namespace GiNaC;

// some very important constants:
const numeric twentyone(21);
const numeric ten(10);
const numeric five(5);

int main()
{
    numeric answer = twentyone;

    answer /= five;
    cout << answer.is_integer() << endl;  // false, it's 21/5
    answer *= ten;
    cout << answer.is_integer() << endl;  // true, it's 42 now!

    return 0;
}