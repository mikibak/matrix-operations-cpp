#pragma clang diagnostic push
#pragma ide diagnostic ignored "modernize-use-auto"
#include <iostream>
#include <cmath>
#include "Matrix.h"
#include "BandMatrix.h"
#include <chrono>
#include <fstream>

#define MY_SIZE 5//968

char filenumber = '0';
char path[] = "/home/mikibak/matrix-operations-cpp/plots/error0.csv";
//47

Matrix LuDecomposition(const Matrix& A, const Matrix& b) {
    int N = A.GetSizeX();

    Matrix L(N, N, 0);
    Matrix U(N, N, 0);

    for (int i = 0; i < N; i++)
        L.SetElement(i, i, 1.0);

    /*for (int j = 0; j < N; j++) {
        for (int i = 0; i <= j; i++) {
            double value = A.GetElement(j, i);
            U.SetElement(j, i, U.GetElement(j, i) + value);
            //U[i][j] += A[i][j];
            for (int k = 0; k <= i - 1; k++) {
                double value1 = L.GetElement(k, i) * U.GetElement(j, k);
                U.SetElement(j, i, U.GetElement(j, i) - value1);
            }
            //U[i][j] -= L[i][k] * U[k][j];

        }

        for (int i = j+1; i < N; i++) {
            double Lij = L.GetElement(j, i);
            for (int k = 0; k <= j - 1; k++) {
                L.SetElement(j, i, Lij - L.GetElement(k, i) * U.GetElement(j, k));
                //L[i][j] -= L[i][k] * U[k][j];
            }
            L.SetElement(j, i, Lij + A.GetElement(j, i));
            L.SetElement(j, i, Lij / U.GetElement(j, j));
            //L[i][j] += A[i][j];
            //L[i][j] /= U[j][j];
        }
    }*/

    //Doolittle
    for (int i = 0; i < N; i++)
    {
        // Upper Triangular
        for (int k = i; k < N; k++)
        {
            // Summation of L(i, j) * U(j, k)
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L.GetElement(j, i) * U.GetElement(k, j));

            // Evaluating U(i, k)
            U.SetElement(k, i, A.GetElement(k, i) - sum);
        }

        // Lower Triangular
        for (int k = i; k < N; k++)
        {
            if (i == k)
                L.SetElement(i, i, 1); // Diagonal as 1
            else
            {
                // Summation of L(k, j) * U(j, i)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L.GetElement(j, k) * U.GetElement(i, j));

                // Evaluating L(k, i)
                L.SetElement(i, k, (A.GetElement(i, k) - sum) / U.GetElement(i, i));
            }
        }
    }

    L.Print();
    U.Print();
    (L*U).Print();

    /* solving equation L * y = b for y using forward substitution method */
    Matrix x = Matrix(1, A.GetSizeY(), 0);
    Matrix y = Matrix(1, A.GetSizeY(), 0);

    for (int i = 0; i < N; i++)
    {
        double S = 0;

        for (int j = 0; j < i; j++) {
            S += L.GetElement(i, j) * y.GetElement(0, j);
        }
        //S += L->A[i][j] * y[j];

        //y[i] = (b[i] - S) / L->A[i][i];
        double value = (b.GetElement(0, 1) - S) / L.GetElement(i, i);
        y.SetElement(0, i, value);
    }

    /* solving equation U * x = y using backward substitution method */
    for (int i = N - 1; i >= 0; i--)
    {
        double S = 0;

        for (int j = i + 1; j < N; j++) {
            S += U.GetElement(i, j);
        } //S += U->A[i][j] * x[j];

        double value = (y.GetElement(0, i) - S) / U.GetElement(i, i);
        x.SetElement(0, i, value);
        //x[i] = (y[i] - S) / U->A[i][i];
    }

    /*  //Podstawienie w przód dla Ly = b
      for (int i = 0; i < N; i++) {
          double val = b.GetElement(i, 0);//b[i][0];
          for (int j = 0; j < i; j++) {
              if (j != i) val -= L.GetElement(i, j) * y.GetElement(j, 0);//L[i][j] * y[j][0];
          }

          y.SetElement(i, 0, val/L.GetElement(i, i));
          //y[i][0] = val / L[i][i];
      }

      //Podstawienie wstecz dla Ux = y
      for (int i = N - 1; i >= 0; i--) {
          double val = y.GetElement(i, 0);//y[i][0];
          for (int j = i; j < N; j++) {
              if (j != i) val -= U.GetElement(i, j) * x.GetElement(j, 0);//U[i][j] * x[j][0];
          }

          x.SetElement(i, 0, val/U.GetElement(i, i));
          //x[i][0] = val / U[i][i];
      }*/

    double error = (A*x - b).Norm();
    std::cout << "Residual error norm for LU decomposition: " << error << std::endl;

    return x;

    //auto end = chrono::high_resolution_clock::now();
    //auto difference = end - start;
    //*duration = chrono::duration<double, milli>(difference).count();

    //double* r = new double[N];		//residual vector
    //residual(A, b, x, r);
    //*n = norm(r);
}


Matrix Jacobi(const Matrix& A, const Matrix& b, int max_iter, double tol, bool save) {

    double error[max_iter];
    for(int i = 0; i < max_iter; i++) {
        error[i] = 0;
    }

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

        error[k] = (A*r - b).Norm();

        //std::cout << error[k];

        if(error[k] < tol) {
            break;
        }

        if(k == max_iter - 1) std::cout << "Reached max iter" << std::endl;
    }

    if(A.GetSizeX() == MY_SIZE) std::cout << "Number of iterations for Jacobi: " << k << std::endl;
    //save error for plot

    if(save) {
        bool comma = false;
        path[47] = filenumber;
        std::ofstream out2(path);
        for (auto col : error) {
            if (comma) out2 << ',';
            comma = true;
            out2 << col;
        }
        filenumber++;
    }

    return r;
}


Matrix GaussSeidel(const Matrix& A, const Matrix& b, int max_iter, double tol, bool save) {

    double error[max_iter];
    for(int i = 0; i < max_iter; i++) {
        error[i] = 0;
    }

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

        error[k] = (A*r - b).Norm();

        if(error[k] < tol) {
            break;
        }

        if(k == max_iter - 1) std::cout << "Reached max iter" << std::endl;
    }

    if(A.GetSizeX() == MY_SIZE) std::cout << "Number of iterations for Gauss-Seidel: " << k << std::endl;

    //save error for plot

    if(save) {
        bool comma = false;
        path[47] = filenumber;
        std::ofstream out2(path);
        for (auto col : error) {
            if (comma) out2 << ',';
            comma = true;
            out2 << col;
        }
        filenumber++;
    }

    return r;
}

long TimeWrapJacobi(const Matrix& A, const Matrix& b, int max_iter, double tol) {
    auto start = std::chrono::high_resolution_clock::now();
    Jacobi(A, b, max_iter, tol, false);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    return duration.count();
}

long TimeWrapGaussSeidel(const Matrix& A, const Matrix& b, int max_iter, double tol) {
    auto start = std::chrono::high_resolution_clock::now();
    GaussSeidel(A, b, max_iter, tol, false);

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    return duration.count();
}

void ValuesForPlots(Matrix* A, Matrix* b, int max_iter, double tol, int a1, int a2, int a3) {
    //N = {100, 500, 1000, 2000, 3000 . . . }

    //int N[] = {10, 50, 100, 200, 300};
    int N[] = {100, 500, 1000, 2000, 3000};

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
        std::cout << i << ": " << GaussSeidelIter[i] << " ms" << std::endl;
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
        std::cout << i << ": " << JacobiIter[i] << " ms" << std::endl;
    }

    //export to csv
    bool comma = false;
    std::ofstream out("/home/mikibak/matrix-operations-cpp/GaussSeidel.csv");
    for (auto col : GaussSeidelIter) {
        if (comma) out << ',';
        comma = true;
        out << col;
    }

    comma = false;
    std::ofstream out2("/home/mikibak/matrix-operations-cpp/Jacobi.csv");
    for (auto col : JacobiIter) {
        if (comma) out2 << ',';
        comma = true;
        out2 << col;
    }
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
    Matrix x = Jacobi(*A, *b, max_iter, 1e-6, true);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << " ms"  << std::endl;
    //x.Print();


    //GAUSS-SEIDEL:

    start = std::chrono::high_resolution_clock::now();
    Matrix x2 = GaussSeidel(*A, *b, max_iter, 1e-6, true);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "Time: " << duration.count() << " ms"  << std::endl;
    //x2.Print();


    ValuesForPlots(A, b, max_iter, 1e-6, a1, a2, a3);


    //C:
    a1 = 3;
    a2 = -1;
    a3 = -1;

    /*A = new BandMatrix(MY_SIZE, MY_SIZE, a1, a2, a3);
    Jacobi(*A, *b, max_iter, 1e-6, true);
    GaussSeidel(*A, *b, max_iter, 1e-6, true);
    std::cout << "Metody nie zbiegają się!" << std::endl;*/

    //LU:

    start = std::chrono::high_resolution_clock::now();
    Matrix LU_result = LuDecomposition(*A, *b);
    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    LU_result.Print();
    std::cout << "LU decomposition time: " << duration.count() << " ms"  << std::endl;

}