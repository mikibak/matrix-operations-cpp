#include <iostream>
#include "Matrix.h"
#include "BandMatrix.h"

#define MY_SIZE 9 * 6 * 8

int main()
{
    std::cout << "Hello World!\n";
    //auto* matrix = new Matrix(MY_SIZE, MY_SIZE);
    Matrix* matrix_ptr = new Matrix(10, 10);
    Matrix matrix = *matrix_ptr;
    matrix.Print();
    //188968
    //Matrix* band_matrix = new BandMatrix(10, 10, 5 + 9, -1, -1);
    Matrix* band_matrix = new BandMatrix(10, 10, 1, 2, 3);
    Matrix* band_matrix2 = new BandMatrix(10, 10, 1, 2, 3);
    band_matrix->Print();
    *band_matrix = *band_matrix + *band_matrix2;
    band_matrix->Print();
}