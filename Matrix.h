//
// Created by mikibak on 4/20/23.
//

#ifndef MATRIX_OPERATIONS_CPP_MATRIX_H
#define MATRIX_OPERATIONS_CPP_MATRIX_H


class Matrix {
protected:
    Matrix();

    int COLS;
    int ROWS;
    int** columns;
    void CheckIfSizeEqual(const Matrix &other);
    void CheckIfValidIndex(int x, int y);


public:
    Matrix(int size_x, int size_y);
    Matrix(int size_x, int size_y, int value);

    int GetElement(int x, int y);
    void SetElement(int x, int y, int value);
    int GetSizeX() const;
    int GetSizeY() const;
    void Print();
    friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
    friend Matrix operator-(const Matrix& lhs, const Matrix& rhs);
    Matrix& operator+=(const Matrix& rhs);

    friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);
    Matrix& operator*=(const Matrix& rhs);

    friend Matrix operator/(const Matrix& m1, const Matrix& m2);

    ~Matrix();

    Matrix *LowerTriangle();

    Matrix *UpperTriangle();

    Matrix *Diagonal();

    double Norm() const;

    Matrix forwardSubstitution(const Matrix& b) const;
};


#endif //MATRIX_OPERATIONS_CPP_MATRIX_H
