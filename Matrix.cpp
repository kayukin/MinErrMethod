//
// Created by kayukin on 03.11.15.
//

#include <string.h>
#include "Matrix.h"

Matrix::Matrix(size_t N) {
    size = N;
    matr = new double *[size];
    for (size_t i = 0; i < size; i++) {
        matr[i] = new double[size];
        memset(matr[i], 0, sizeof(double) * size);
    }
}

Matrix::Matrix(const Matrix &matrix) {
    Copy(matrix);
}

Matrix::~Matrix() {
    Destroy();
}

Matrix &Matrix::operator=(const Matrix &matrix) {
    if (this == &matrix) {
        return *this;
    }
    Destroy();
    Copy(matrix);
}

void Matrix::Destroy() {
    for (size_t i = 0; i < size; i++)
        delete[] matr[i];
    delete[] matr;
}

void Matrix::Copy(const Matrix &matrix) {
    size = matrix.size;
    matr = new double *[size];
    for (size_t i = 0; i < size; i++) {
        matr[i] = new double[size];
        for (size_t j = 0; j < size; j++) {
            matr[i][j] = matrix.matr[i][j];
        }
    }
}

Matrix Matrix::operator*(const Matrix &B) const {
    Matrix result(size);
    for (size_t i = 0; i < size; i++)
        for (size_t j = 0; j < size; j++) {
            for (size_t k = 0; k < size; k++)
                result.matr[i][j] += this->matr[i][k] * B.matr[k][j];
        }
}