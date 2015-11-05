#ifndef MINERRMETHOD_MINIMALERRORSMETHOD_H
#define MINERRMETHOD_MINIMALERRORSMETHOD_H

#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include "Generator.h"

class MinimalErrorsMethod {
private:
public:
    MinimalErrorsMethod();
    static int Solve(const Matrix &, Vector &, Vector &, double);
};

double ABS_ERROR(Vector vector1, Vector vector2);

double ABS_NEV(Matrix p, Vector v, Vector f);

double OTN_ERROR(Vector v1, Vector v2);

double OTN_NEV(Matrix p, Vector v, Vector f);

void Mult_On_Inv_Matr_For_Sopr(Matrix& p, const Matrix& B);

#endif //MINERRMETHOD_MINIMALERRORSMETHOD_H
