#pragma clang diagnostic push
#pragma ide diagnostic ignored "modernize-use-auto"
#include <iostream>
#include <cmath>
#include "Matrix.h"
#include "BandMatrix.h"

#define MY_SIZE 5//(9 * 6 * 8)

Matrix Jacobi(const Matrix& A, const Matrix& b, const Matrix& L, const Matrix& U, const Matrix& D, int max_iter, double tol) {
    Matrix x = Matrix(1, MY_SIZE, 1);
    // Jacobi iteration
    for (int k = 0; k < max_iter; k++) {
        // Compute the next approximation
        Matrix minus_D = Matrix(D.GetSizeX(), D.GetSizeY(), -1) * D;
        Matrix LU = L + U;
        LU.Print();
        x.Print();
        Matrix LUx = LU * x;
        LUx.Print();
        Matrix sub1 = minus_D.forwardSubstitution(LUx);
        sub1.Print();
        Matrix sub2 = D.forwardSubstitution(b);
        sub2.Print();
        Matrix x_next = sub1 + sub2;

        // Check if the difference is smaller than the tolerance
        Matrix m = x_next - x;
        if (m.Norm() < tol) {
            return x_next;
        }

        // Update x
        x = x_next;
    }

    // Return the last approximation
    return x;
}


int main()
{
    std::cout << "Hello World!\n";
    //auto* matrix = new Matrix(MY_SIZE, MY_SIZE);
    //Matrix* matrix_ptr = new Matrix(10, 10);
    //Matrix matrix = *matrix_ptr;
    //matrix.Print();
    //188968
    //Matrix* band_matrix = new BandMatrix(10, 10, 5 + 9, -1, -1);
/*
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
*/


    //TASK A

    int a1 = 5 + 9;
    int a2 = -1;
    int a3 = -1;

    int max_iter = 1000;

    Matrix* A = new BandMatrix(MY_SIZE, MY_SIZE, a1, a2, a3);
    Matrix* b = new Matrix(1, MY_SIZE);

    for(int i = 0; i < MY_SIZE; i++) {
        int value = sin( i * (8 + 1));
        b->SetElement(0, i, value);
    }

    Matrix* L = A->LowerTriangle();
    Matrix* U = A->UpperTriangle();
    Matrix* D = A->Diagonal();

    Matrix x = Jacobi(*A, *b, *L, *U, *D, 1000, 1e-6);

    //while(true) {

        //r Matrix x_next = -D\(L + U) * r + D\b;
        //error = norm(M*r - b)

        //if(error < pow(10, -14)) {
        //    break;
        //}
    //}
}