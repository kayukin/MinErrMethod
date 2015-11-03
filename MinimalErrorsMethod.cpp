#include <cmath>
#include "MinimalErrorsMethod.h"
#include "Generator.h"

MinimalErrorsMethod::MinimalErrorsMethod() {

}

int MinimalErrorsMethod::Solve(const Matrix &A, Vector &xk, Vector &f, double eps) {
    int numberOfIterations = 0;
    Vector rk(A.getSize());
    Vector ark(A.getSize());
    double arkrk = 0;
    double rkrk = 0;
    double tk = 0;

    do {
        //r(k)=A*x(k)-f

        rk = A * xk;
        rk = rk - f;

        //tk=(r(k),r(k))/(A*r(k),r(k));

        arkrk = 0;
        ark = A * rk;
        for (int i = 0; i < A.getSize(); i++)
            arkrk += ark(i) * rk(i);
        rkrk = 0;
        for (int i = 0; i < A.getSize(); i++)
            rkrk += rk(i) * rk(i);
        tk = rkrk / arkrk;

        // x(k+1)=x(k)-t(k)*r(k),

        for (int i = 0; i < A.getSize(); i++)
            xk(i) -= tk * rk(i);

        numberOfIterations++;
    }
    while (rk.Norma() > 1e-6 && rk.Norma() < 1e+30 && numberOfIterations < (10e+7) / A.getSize());
    return numberOfIterations;
}

void MinimalErrorsMethod::main() {
    double min = 1;
    double max = 0;
    int n;
    int iter;
    cout << "‚ўҐ¤ЁвҐ а §¬Ґа ¬ ваЁжл -> ";
    //cin >> n;
    n = 10;
    Generator generator(n, min, max);
    Vector x(n);//вектор точного решения
    //вектор точных значений
    Matrix p = generator.getA();//исходная матрица
    Matrix B = p;//исходная матрица
    //vector(n, x);//точный вектор
    //Jordan_Form(n, p, min, max);
    //cout << fabs(Obuslov(max, min)) << endl;
    double eps = 1e-5;
    //cin >> eps;
    B = p.Sopr();
    p = B * p;//Matr_On_Matr2(n, B, p);
    //Mult_On_Inv_Matr_For_Sopr(n, p, B);
    Vector b = p * x;//вектор значений
    Vector a(n);//приближенное решение
    iter = Solve(p, a, b, eps);
    cout << "===========================================" << endl;
    cout << "||  ЂЎб.®иЁЎЄ    ||" << ABS_ERROR(a, x) << "||" << endl;
    cout << "===========================================" << endl;
    cout << "===========================================" << endl;
    cout << "||  ЂЎб.­Ґўп§Є   ||" << ABS_NEV(p, a, b) << "||" << endl;
    cout << "===========================================" << endl;
    cout << "===========================================" << endl;
    cout << "||  Ћв­.®иЁЎЄ    ||" << OTN_ERROR(a, x) << "||" << endl;
    cout << "===========================================" << endl;
    cout << "===========================================" << endl;
    cout << "||  Ћв­.­Ґўп§Є   ||" << OTN_NEV(p, a, b) << "||" << endl;
    cout << "===========================================" << endl;
}
