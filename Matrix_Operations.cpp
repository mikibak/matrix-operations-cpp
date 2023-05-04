#pragma clang diagnostic push
#pragma ide diagnostic ignored "modernize-use-auto"
#include <iostream>
#include <cmath>
#include "Matrix.h"
#include "BandMatrix.h"
#include <chrono>

#define MY_SIZE 968



Matrix Jacobi(const Matrix& A, const Matrix& b, int max_iter, double tol) {

    Matrix L = A.LowerTriangle();
    Matrix U = A.UpperTriangle();
    Matrix D = A.Diagonal();

    Matrix r = Matrix(1, A.GetSizeX(), 1);
    Matrix minus_D = Matrix(D.GetSizeX(), D.GetSizeY());
    Matrix LU = L + U;
    Matrix DforwardSubstitutionb =  D.forwardSubstitution(b);

    for(int i = 0; i < A.GetSizeX(); i++) {
        minus_D.SetElement(i, i, -1 *D.GetElement(i, i));
    }
    // Jacobi iteration
    int k;
    for (k = 0; k < max_iter; k++) {

        r = minus_D.forwardSubstitution(LU * r) + DforwardSubstitutionb;

        double error = (A*r - b).Norm();

        if(error < tol) {
            break;
        }
    }

    if(A.GetSizeX() == MY_SIZE) std::cout << "Number of iterations for Jacobi: " << k << std::endl;
    return r;
}


Matrix GaussSeidel(const Matrix& A, const Matrix& b, int max_iter, double tol) {

    Matrix L = A.LowerTriangle();
    Matrix U = A.UpperTriangle();
    Matrix D = A.Diagonal();

    Matrix r = Matrix(1, A.GetSizeX(), 1);
    Matrix minus_DL = Matrix(D.GetSizeX(), D.GetSizeY());
    Matrix DLforwardSubstitutionb =  (D + L).forwardSubstitution(b);

    //create DL
    for(int i = 0; i < A.GetSizeX(); i++) {
        for(int j = 0; j < A.GetSizeX(); j++) {
            double value = -1 * (D.GetElement(i, j) + L.GetElement(i, j));
            minus_DL.SetElement(i, j, value);
        }
    }
    // Jacobi iteration
    int k;
    for (k = 0; k < max_iter; k++) {

        //r = -(D+L) \ (U*r) + (D+L) \ b;
        r = minus_DL.forwardSubstitution(U * r) + DLforwardSubstitutionb;

        double error = (A*r - b).Norm();

        if(error < tol) {
            break;
        }

    }

    if(A.GetSizeX() == MY_SIZE) std::cout << "Number of iterations for Gauss-Seidel: " << k << std::endl;
    return r;
}

long TimeWrapJacobi(const Matrix& A, const Matrix& b, int max_iter, double tol) {
    auto start = std::chrono::high_resolution_clock::now();
    Jacobi(A, b, max_iter, tol);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    return duration.count();
}

long TimeWrapGaussSeidel(const Matrix& A, const Matrix& b, int max_iter, double tol) {
    auto start = std::chrono::high_resolution_clock::now();
    GaussSeidel(A, b, max_iter, tol);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    return duration.count();
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

    //JACOBI:

    auto start = std::chrono::high_resolution_clock::now();

    Matrix x = Jacobi(*A, *b, max_iter, 1e-6);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << " ms"  << std::endl;
    //x.Print();


    //GAUSS-SEIDEL:

    start = std::chrono::high_resolution_clock::now();

    Matrix x2 = GaussSeidel(*A, *b, max_iter, 1e-6);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << " ms"  << std::endl;
    //x2.Print();


    //C:
    a1 = 3;
    a2 = -1;
    a3 = -1;

    A = new BandMatrix(MY_SIZE, MY_SIZE, a1, a2, a3);
    Jacobi(*A, *b, max_iter, 1e-6);
    GaussSeidel(*A, *b, max_iter, 1e-6);

    //N = {100, 500, 1000, 2000, 3000 . . . }
    std::cout << "\n\n";
    A = new BandMatrix(100, 100, a1, a2, a3);
    b = new Matrix(1, 100);
    for(double i = 0; i < 100; i++) {
        double value = sin( i * (8 + 1));
        b->SetElement(0, i, value);
    }
    long GaussSeidel100 = TimeWrapGaussSeidel(*A, *b, max_iter, 1e-6);
    std::cout << "100: " << GaussSeidel100 << " s" << std::endl;

    A = new BandMatrix(500, 500, a1, a2, a3);
    b = new Matrix(1, 500);
    for(double i = 0; i < 500; i++) {
        double value = sin( i * (8 + 1));
        b->SetElement(0, i, value);
    }
    long GaussSeidel500 = TimeWrapGaussSeidel(*A, *b, max_iter, 1e-6);
    std::cout << "500: " << GaussSeidel500 << " s" << std::endl;

    A = new BandMatrix(1000, 1000, a1, a2, a3);
    b = new Matrix(1, 1000);
    for(double i = 0; i < 1000; i++) {
        double value = sin( i * (8 + 1));
        b->SetElement(0, i, value);
    }
    long GaussSeidel1000 = TimeWrapGaussSeidel(*A, *b, max_iter, 1e-6);
    std::cout << "1000: " << GaussSeidel1000 << " s" << std::endl;

    for(int i = 1000; i < 4000; i+= 1000) {
        A = new BandMatrix(i, i, a1, a2, a3);
        b = new Matrix(1, i);
        for(double j = 0; j < j; j++) {
            double value = sin( j * (8 + 1));
            b->SetElement(0, j, value);
        }
        long GaussSeidel500 = TimeWrapGaussSeidel(*A, *b, max_iter, 1e-6);
        std::cout << i << ": " << GaussSeidel500 << " s" << std::endl;
    }


    std::cout << "\n\n";
    A = new BandMatrix(100, 100, a1, a2, a3);
    b = new Matrix(1, 100);
    for(double i = 0; i < 100; i++) {
        double value = sin( i * (8 + 1));
        b->SetElement(0, i, value);
    }
    long Jacobi100 = TimeWrapJacobi(*A, *b, max_iter, 1e-6);
    std::cout << "100: " << Jacobi100 << " s" << std::endl;

    A = new BandMatrix(500, 500, a1, a2, a3);
    b = new Matrix(1, 500);
    for(double i = 0; i < 500; i++) {
        double value = sin( i * (8 + 1));
        b->SetElement(0, i, value);
    }
    long Jacobi500 = TimeWrapJacobi(*A, *b, max_iter, 1e-6);
    std::cout << "500: " << Jacobi500 << " s" << std::endl;

    for(int i = 1000; i < 4000; i+= 1000) {
        A = new BandMatrix(i, i, a1, a2, a3);
        b = new Matrix(1, i);
        for(double i = 0; i < i; i++) {
            double value = sin( i * (8 + 1));
            b->SetElement(0, i, value);
        }
        long Jacobi500 = TimeWrapJacobi(*A, *b, max_iter, 1e-6);
        std::cout << i << ": " << Jacobi500 << " s" << std::endl;
    }
}