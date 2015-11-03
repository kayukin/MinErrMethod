#include <iostream>
#include "Matrix.h"

using namespace std;

int main() {
    Matrix A(3);

    int k = 1;
    for (int i = 0; i < A.getSize(); i++) {
        for (int j = 0; j < A.getSize(); ++j) {
            A(i, j) = k++;
        }
    }
    A(0, 0) = 10;
    cout << A << endl;
    cout << A.Determinant();
    return 0;
}