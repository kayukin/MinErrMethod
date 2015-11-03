//
// Created by kayukin on 03.11.15.
//

#ifndef MINERRMETHOD_MATRIX_H
#define MINERRMETHOD_MATRIX_H


#include <stddef.h>

class Matrix {
protected:
    double **matr;
    size_t size;

    void Destroy();

    void Copy(const Matrix &);

public:
    Matrix(size_t N);

    Matrix(const Matrix &);

    ~Matrix();

    Matrix &operator=(const Matrix &);

    Matrix operator*(const Matrix &B) const;

};


#endif //MINERRMETHOD_MATRIX_H
