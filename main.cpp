#include <iostream>
#include "MinimalErrorsMethod.h"

using namespace std;

int main() {
    double beta = 10;
    double alpha = 1;
    size_t n = 10;
    double eps = 1e-5;
    int iter;
    Generator generator(n, alpha, beta);
    Vector x(n);
    Matrix p = generator.getA();
    for (size_t i = 0; i < x.getSize(); i++)
        x(i) = 1;
    Vector b = p * x;
    Vector a(n);
    iter = MinimalErrorsMethod::Solve(p, a, b, eps);
    cout << "================================" << endl;
    cout << "||  ABS_ERROR || = " << ABS_ERROR(a, x) << endl;
    cout << "================================" << endl;
    cout << "================================" << endl;
    cout << "||  ABS_NEV   || = " << ABS_NEV(p, a, b) << endl;
    cout << "================================" << endl;
    cout << "================================" << endl;
    cout << "||  OTN_ERROR || = " << OTN_ERROR(a, x) << endl;
    cout << "================================" << endl;
    cout << "================================" << endl;
    cout << "||  OTN_NEV   || = " << OTN_NEV(p, a, b) << endl;
    cout << "================================" << endl;
    cout << "================================" << endl;
    cout << "||  OBUSLOVLENNOST'   || = " << generator.getObusl() << endl;
    cout << "================================" << endl;
}