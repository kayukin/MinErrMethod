#include <cmath>
#include "Vector.h"

Vector::Vector(size_t N) {
    size = N;
    vector = new double[size];
    memset(vector, 0, sizeof(double) * size);
}

Vector::Vector(const Vector &vector) {
    this->size = vector.size;
    this->vector = new double[size];
    for (int i = 0; i < this->size; ++i) {
        this->vector[i] = vector.vector[i];
    }
}

Vector &Vector::operator=(const Vector &vector) {
    if (this == &vector)
        return *this;
    delete[] this->vector;
    this->size = vector.size;
    this->vector = new double[size];
    for (int i = 0; i < this->size; ++i) {
        this->vector[i] = vector.vector[i];
    }
    return *this;
}

Vector::~Vector() {
    delete[] vector;
}

ostream &operator<<(ostream &os, const Vector &vector) {
    for (int i = 0; i < vector.size; ++i) {
        os << vector.vector[i] << ' ';
    }
    return os;
}

double &Vector::operator()(size_t i) {
    return vector[i];
}

double Vector::operator()(size_t i) const {
    return vector[i];
}

Vector Vector::operator-(const Vector &vector) const {
    Vector result(*this);
    for (size_t i = 0; i < result.getSize(); i++) {
        result(i) -= vector(i);
    }
    return result;
}

size_t Vector::getSize() const {
    return size;
}

double Vector::Norma() const {
    double norma = 0;
    for (int i = 0; i < size; i++)
        norma += vector[i] * vector[i];
    return sqrt(norma);
}
