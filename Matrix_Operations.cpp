#include <iostream>
#include <math.h>
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


    //TASK A

    int a1 = a1 = 5 + 9;
    int a2 = -1;
    int a3 = -1;

    Matrix A = BandMatrix(MY_SIZE, MY_SIZE, a1, a2, a3);

    Matrix b = Matrix(1, MY_SIZE);
    for(int i = 0; i < MY_SIZE; i++) {
        b.SetElement(0, i, sin( i * (8 + 1)));
    }

}