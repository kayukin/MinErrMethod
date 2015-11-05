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
    return *this;
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
    return result;
}

double Matrix::Norm() const {
    return 0;
}

double Matrix::Determinant() const {
    int l;
    double d;
    double sum11 = 1, sum12 = 0, sum21 = 1, sum22 = 0;
    for (int i = 0; i < size; i++) {
        sum11 = 1;
        l = (int) (2 * size - 1 - i);
        sum21 = 1;
        for (int j = 0; j < size; j++) {
            sum21 *= matr[j][l % size];
            l--;
            sum11 *= matr[j][(j + i) % (size)];
        }
        sum22 += sum21;
        sum12 += sum11;
    }
    d = sum12 - sum22;
    return d;
}

std::ostream &operator<<(std::ostream &os, const Matrix &matrix) {
    for (size_t i = 0; i < matrix.size; i++) {
        for (size_t j = 0; j < matrix.size; j++)
            os << matrix.matr[i][j] << ' ';
        os << std::endl;
    }
    return os;
}

double &Matrix::operator()(size_t i, size_t j) {
    return matr[i][j];
}

double Matrix::operator()(size_t i, size_t j) const {
    return matr[i][j];
}

size_t Matrix::getSize() const {
    return size;
}

Matrix Matrix::Sopr() {
    Matrix B(size);
    for (int i = 0; i < size - 1; i++)
        B(0, i) = 1 / sqrt((i + 1) * (i + 2));
    //B(0, i) = 1 / sqrt(size);
    for (int i = 1; i < size - 1; i++) {
        B(i, i - 1) = -sqrt(i) / sqrt(i + 1);
        for (int j = i; j < size - 1; j++)
            B(i, j) = 1 / sqrt((j + 1) * (j + 2));
        //B(i, j) = 1 / sqrt(size);
    }
    B(size - 1, size - 2) = -1 / sqrt(size * (size - 1));
    B(size - 1, size - 1) = 1 / sqrt(size);
    return B;
}