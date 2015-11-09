#include <iostream>
#include <fstream>
#include "MinimalErrorsMethod.h"

using namespace std;

void GenerateAndSolve(double alpha, double beta, size_t n, double eps, ostream &os);

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
    Vector b = p * x;
    Vector a(n);
    MinimalErrorsMethod::Solve(p, a, b, eps);
    cout << "|| A || = " << generator.getNorm();
    cout << "|| A_INV || = " << generator.getNorm_inv();

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
    /*ofstream file("output.txt");
    file << "alpha  " << "beta   " << "||A||       " << "||A_INV||      " << "ABS_ERROR       " << "ABS_NEV     " <<
    "OTN_ERROR      " << "OTN_NEV  " << "OBUSLOVLENNOST' " << endl;
    for (beta = 10; beta <= 1e16; beta *= beta) {
        GenerateAndSolve(alpha, beta, n, eps, file);
    }
    beta=1;
    for(alpha=1;alpha>=1e-16;alpha/=alpha){
        GenerateAndSolve(alpha, beta, n, eps, file);*/

}

void GenerateAndSolve(double alpha, double beta, size_t n, double eps, ostream &os) {
    Generator generator(n, alpha, beta);
    Vector x(n);
    Matrix p = generator.getA();
    for (size_t i = 0; i < x.getSize(); i++)
        x(i) = 1;
    Vector b = p * x;
    Vector a(n);
    MinimalErrorsMethod::Solve(p, a, b, eps);
    os << alpha << "|" << beta << "|" << generator.getNorm() << "   |     " << generator.getNorm_inv() << "  |  " <<
    ABS_ERROR(a, x) << "|" << ABS_NEV(p, a, b) << "|" << OTN_ERROR(a, x) << "|" << OTN_NEV(p, a, b) << "|" <<
    generator.getObusl() << endl;
}