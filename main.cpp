#include <iostream>
#include "Matrix.h"
#include "Generator.h"

using namespace std;

int main() {
    Generator generator(10, 1, 5);
    cout << generator.getA();
    return 0;
}