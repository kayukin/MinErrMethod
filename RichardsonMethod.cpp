#include "RichardsonMethod.h"

Vector RichardsonMethod::Solve(Matrix A, Vector f, double T) {
    Vector U(A.getSize());
    for (size_t i = 0; i < A.getSize(); i++)
        U(i) = 1;
    while (a > 0) {
        U = U - T * (A * U - f);
    }
    return U;
}

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

void Mult_On_Inv_Matr_For_Sopr(Matrix &p, const Matrix &B) {

}
