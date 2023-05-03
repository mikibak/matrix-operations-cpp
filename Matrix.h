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
    double** columns;
    void CheckIfSizeEqual(const Matrix &other) const;
    void CheckIfValidIndex(int x, int y) const;


public:
    Matrix(int size_x, int size_y);
    Matrix(int size_x, int size_y, int value);

    double GetElement(int x, int y) const;
    void SetElement(int x, int y, double value);
    int GetSizeX() const;
    int GetSizeY() const;
    void Print() const;
    friend Matrix operator+(const Matrix& lhs, const Matrix& rhs);
    friend Matrix operator-(const Matrix& lhs, const Matrix& rhs);
    Matrix& operator+=(const Matrix& rhs);

    friend Matrix operator*(const Matrix& lhs, const Matrix& rhs);
    Matrix& operator*=(const Matrix& rhs);

    friend Matrix operator/(const Matrix& m1, const Matrix& m2);

    Matrix& operator=(const Matrix& rhs) noexcept;

    ~Matrix();

    Matrix *LowerTriangle();

    Matrix *UpperTriangle();

    Matrix *Diagonal();

    double Norm() const;

    Matrix forwardSubstitution(const Matrix& b) const;
};


#endif //MATRIX_OPERATIONS_CPP_MATRIX_H
