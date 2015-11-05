#ifndef MINERRMETHOD_MINIMALERRORSMETHOD_H
#define MINERRMETHOD_MINIMALERRORSMETHOD_H

#include "Matrix.h"
#include "Vector.h"
#include <cmath>
#include "Generator.h"

class MinimalErrorsMethod {
private:

    static int Solve(const Matrix &, Vector &, Vector &, double);

public:
    MinimalErrorsMethod();

    static void main();

};

double ABS_ERROR(Vector vector1, Vector vector2);

double ABS_NEV(Matrix p, Vector v, Vector f);

double OTN_ERROR(Vector v1, Vector v2);

double OTN_NEV(Matrix p, Vector v, Vector f);

#endif //MINERRMETHOD_MINIMALERRORSMETHOD_H
