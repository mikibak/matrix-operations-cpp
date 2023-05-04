#pragma clang diagnostic push
#pragma ide diagnostic ignored "modernize-use-auto"
#include <iostream>
#include <cmath>
#include "Matrix.h"
#include "BandMatrix.h"
#include <chrono>
#include <fstream>

#define MY_SIZE 968

Matrix LuDecomposition(const Matrix& A, const Matrix& b) {
    Matrix U = A.UpperTriangle() + A.Diagonal();
    Matrix L = A.LowerTriangle() + A.Diagonal();

    int N = A.GetSizeX();

    //auto start = chrono::high_resolution_clock::now();

    /* creating matrices L and U, such that A = L * U */
    for (int i = 0; i < N - 1; ++i)		//columns
    {
        for (int j = i + 1; j < N; ++j)	//rows
        {
            //L->A[j][i] = U->A[j][i] / U->A[i][i];
            double value = U.GetElement(j, i) / U.GetElement(i, i);
            L.SetElement(j, i, value);

            for (int k = i; k < N; ++k) {
                //U->A[j][k] = U->A[j][k] - L->A[j][i] * U->A[i][k];
                double v = U.GetElement(j, k) - L.GetElement(j, i) * U.GetElement(i, k);
                U.SetElement(j, k, v);
            }
        }
    }
    /* now the main equation to solve may be defined as L * U * x = b  */

    /* solving equation L * y = b for y using forward substitution method */
    Matrix x = Matrix(1, A.GetSizeY(), 0);
    Matrix y = Matrix(1, A.GetSizeY(), 0);

    for (int i = 0; i < N; ++i)
    {
        double S = 0;

        for (int j = 0; j < i; ++j) {
            S += L.GetElement(i, j) * y.GetElement(0, j);
        }
            //S += L->A[i][j] * y[j];

        //y[i] = (b[i] - S) / L->A[i][i];
        double value = (b.GetElement(0, 1) - S) / L.GetElement(i, i);
        y.SetElement(0, i, value);
    }

    /* solving equation U * x = y using backward substitution method */
    for (int i = N - 1; i >= 0; --i)
    {
        double S = 0;

        for (int j = i + 1; j < N; ++j) {
            S += U.GetElement(i, j);
        } //S += U->A[i][j] * x[j];

        double value = (y.GetElement(0, i) - S) / U.GetElement(i, i);
        x.SetElement(0, i, value);
        //x[i] = (y[i] - S) / U->A[i][i];
    }

    return x;

    //auto end = chrono::high_resolution_clock::now();
    //auto difference = end - start;
    //*duration = chrono::duration<double, milli>(difference).count();

    //double* r = new double[N];		//residual vector
    //residual(A, b, x, r);
    //*n = norm(r);
}

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
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    return duration.count();
}

long TimeWrapGaussSeidel(const Matrix& A, const Matrix& b, int max_iter, double tol) {
    auto start = std::chrono::high_resolution_clock::now();
    GaussSeidel(A, b, max_iter, tol);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    return duration.count();
}

void ValuesForPlots(Matrix* A, Matrix* b, int max_iter, double tol, int a1, int a2, int a3) {
    //N = {100, 500, 1000, 2000, 3000 . . . }

    int N[] = {10, 50, 100, 200, 300};
    //int N[] = {100, 500, 1000, 2000, 3000};

    long GaussSeidelIter[5];
    long JacobiIter[5];

    std::cout << "\n\n";

    for(int i = 0; i < 5; i++) {
        A = new BandMatrix(N[i], N[i], a1, a2, a3);
        b = new Matrix(1, N[i]);
        for(double j = 0; j < N[i]; j++) {
            double value = sin( j * (8 + 1));
            b->SetElement(0, j, value);
        }
        GaussSeidelIter[i] = TimeWrapGaussSeidel(*A, *b, max_iter, 1e-6);
        std::cout << i << ": " << GaussSeidelIter[i] << " s" << std::endl;
    }

    std::cout << "\n";

    for(int i = 0; i < 5; i++) {
        A = new BandMatrix(N[i], N[i], a1, a2, a3);
        b = new Matrix(1, N[i]);
        for(double j = 0; j < N[i]; j++) {
            double value = sin( j * (8 + 1));
            b->SetElement(0, j, value);
        }
        JacobiIter[i] = TimeWrapJacobi(*A, *b, max_iter, 1e-6);
        std::cout << i << ": " << JacobiIter[i] << " s" << std::endl;
    }

    //export to csv
    std::ofstream out("/home/mikibak/matrix-operations-cpp/GaussSeidel.csv");
    for (auto col : GaussSeidelIter)
        out << col <<',';
    out << '\n';

    std::ofstream out2("/home/mikibak/matrix-operations-cpp/Jacobi.csv");
    for (auto col : JacobiIter)
        out2 << col <<',';
    out2 << '\n';
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
    /*a1 = 3;
    a2 = -1;
    a3 = -1;

    A = new BandMatrix(MY_SIZE, MY_SIZE, a1, a2, a3);
    Jacobi(*A, *b, max_iter, 1e-6);
    GaussSeidel(*A, *b, max_iter, 1e-6);*/

    ValuesForPlots(A, b, max_iter, 1e-6, a1, a2, a3);

}