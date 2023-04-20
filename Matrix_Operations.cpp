#include <iostream>
#include "Matrix.h"
#include "BandMatrix.h"

#define MY_SIZE 9 * 6 * 8

int main()
{
    std::cout << "Hello World!\n";
    //auto* matrix = new Matrix(MY_SIZE, MY_SIZE);
    //Matrix* matrix_ptr = new Matrix(10, 10);
    //Matrix matrix = *matrix_ptr;
    //matrix.Print();
    //188968
    //Matrix* band_matrix = new BandMatrix(10, 10, 5 + 9, -1, -1);
    Matrix* band_matrix = new BandMatrix(10, 10, 1, 2, 3);
    Matrix* band_matrix2 = new BandMatrix(10, 10, 1, 2, 3);
    Matrix matrix = *band_matrix;
    Matrix matrix2 = *band_matrix2;
    matrix.Print();
    matrix += matrix2;
    matrix.Print();

    Matrix m1 = Matrix(10, 10);
    Matrix m2 = Matrix(10, 11);
    m2.Print();
    m1 += m2;
    m1.Print();
}