//
// Created by mikibak on 4/20/23.
//

#ifndef MATRIX_OPERATIONS_CPP_MATRIX_H
#define MATRIX_OPERATIONS_CPP_MATRIX_H


class Matrix {
protected:
    int SIZE_X;
    int SIZE_Y;
    int** columns;


public:
    Matrix(int size_x, int size_y);
    void Print();
    friend Matrix operator+(Matrix lhs, const Matrix& rhs);
    Matrix& operator+=(const Matrix& rhs);

    ~Matrix();

    void CheckIfSizeEqual(const Matrix &other);
};


#endif //MATRIX_OPERATIONS_CPP_MATRIX_H
