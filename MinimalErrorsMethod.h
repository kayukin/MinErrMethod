#ifndef MINERRMETHOD_MINIMALERRORSMETHOD_H
#define MINERRMETHOD_MINIMALERRORSMETHOD_H

#include "Matrix.h"
#include "Vector.h"

class MinimalErrorsMethod {
private:

    static int Solve(const Matrix &, Vector &, Vector &, double);

public:
    MinimalErrorsMethod();

    static void main();

};

double ABS_ERROR(Vector vector1, Vector vector2) {
    Vector v = vector1 - vector2;
    double x = v.Norma();
    return x;
}

double ABS_NEV(Matrix p, Vector v, Vector f) {
    Vector v3 = p * v;
    Vector v4 = v3 - f;
    double x = v4.Norma();
    return x;
}

double OTN_ERROR(Vector v1, Vector v2) {
    double x = ABS_ERROR(v1, v2);
    return x / v2.Norma();
}

double OTN_NEV(Matrix p, Vector v, Vector f) {
    double x = ABS_NEV(p, v, f);
    return x / f.Norma();
}

#endif //MINERRMETHOD_MINIMALERRORSMETHOD_H
