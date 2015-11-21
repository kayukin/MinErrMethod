#include <iostream>
#include "RichardsonMethod.h"

using namespace std;

int main() {
    double beta = 10;
    double alpha = 1;
    size_t n = 10;
    double eps = 1e-5;


    Generator generator(n, alpha, beta);
    Vector x(n);
    Matrix p = generator.getA();
    for (size_t i = 0; i < x.getSize(); i++)
        x(i) = 1;
    Vector b = x;
    double T = 0.0001;
    Vector a = RichardsonMethod::Solve(p, b, T);
    cout << "|| A || = " << generator.getNorm() << endl;
    cout << "|| A_INV || = " << generator.getNorm_inv() << endl;

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