#pragma clang diagnostic push
#pragma ide diagnostic ignored "modernize-use-auto"
#include <iostream>
#include <cmath>
#include "Matrix.h"
#include "BandMatrix.h"

#define MY_SIZE 5//968


Matrix Jacobi(const Matrix& A, const Matrix& b, const Matrix& L, const Matrix& U, const Matrix& D, int max_iter, double tol) {
    Matrix r = Matrix(1, MY_SIZE, 1);
    Matrix minus_D = Matrix(D.GetSizeX(), D.GetSizeY());
    Matrix LU = L + U;
    Matrix DforwardSubstitutionb =  D.forwardSubstitution(b);

    for(int i = 0; i < MY_SIZE; i++) {
        minus_D.SetElement(i, i, -1 *D.GetElement(i, i));
    }
    // Jacobi iteration
    for (int k = 0; k < max_iter; k++) {

        r = minus_D.forwardSubstitution(LU * r) + DforwardSubstitutionb;

        double error = (A*r - b).Norm();

        if(error < tol) {
            std::cout << "Number of iterations for Jacobi: " << k << std::endl;
            return r;
        }
    }

    // Return the last approximation
    return r;
}


Matrix GaussSeidel(const Matrix& A, const Matrix& b, const Matrix& L, const Matrix& U, const Matrix& D, int max_iter, double tol) {
    Matrix r = Matrix(1, MY_SIZE, 1);
    Matrix minus_DL = Matrix(D.GetSizeX(), D.GetSizeY());
    Matrix DLforwardSubstitutionb =  (D + L).forwardSubstitution(b);

    //create DL
    for(int i = 0; i < MY_SIZE; i++) {
        for(int j = 0; j < MY_SIZE; j++) {
            double value = -1 * (D.GetElement(i, j) + L.GetElement(i, j));
            minus_DL.SetElement(i, j, value);
        }
    }
    // Jacobi iteration
    for (int k = 0; k < max_iter; k++) {

        //r = -(D+L) \ (U*r) + (D+L) \ b;
        r = minus_DL.forwardSubstitution(U * r) + DLforwardSubstitutionb;

        double error = (A*r - b).Norm();

        if(error < tol) {
            std::cout << "Number of iterations for Gauss-Seidel: " << k << std::endl;
            return r;
        }

    }

    // Return the last approximation
    return r;
}


int main()
{
    //TASK A

    int a1 = 5 + 9;
    int a2 = -1;
    int a3 = -1;

    int max_iter = 1000;

    Matrix* A = new BandMatrix(MY_SIZE, MY_SIZE, a1, a2, a3);
    Matrix* b = new Matrix(1, MY_SIZE);

    for(double i = 0; i < MY_SIZE; i++) {
        double value = sin( i * (8 + 1));
        b->SetElement(0, i, value);
    }

    Matrix* L = A->LowerTriangle();
    Matrix* U = A->UpperTriangle();
    Matrix* D = A->Diagonal();

    Matrix x = Jacobi(*A, *b, *L, *U, *D, max_iter, 1e-6);
    x.Print();

    Matrix x2 = GaussSeidel(*A, *b, *L, *U, *D, max_iter, 1e-6);
    x2.Print();
    //while(true) {

        //r Matrix x_next = -D\(L + U) * r + D\b;
        //error = norm(M*r - b)

        //if(error < pow(10, -14)) {
        //    break;
        //}
    //}
}