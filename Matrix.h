//
// Created by mikibak on 4/20/23.
//

#ifndef MATRIX_OPERATIONS_CPP_MATRIX_H
#define MATRIX_OPERATIONS_CPP_MATRIX_H


class Matrix {
protected:
    Matrix();

    int SIZE_X;
    int SIZE_Y;
    int** columns;
    void CheckIfSizeEqual(const Matrix &other);
    void CheckIfValidIndex(int x, int y);


public:
    Matrix(int size_x, int size_y);
    int GetElement(int x, int y);
    void SetElement(int x, int y, int value);
    void Print();
    friend Matrix operator+(Matrix lhs, const Matrix& rhs);
    Matrix& operator+=(const Matrix& rhs);

    ~Matrix();

    Matrix *LowerTriangle();

    Matrix *UpperTriangle();

    Matrix *Diagonal();
};


#endif //MATRIX_OPERATIONS_CPP_MATRIX_H
