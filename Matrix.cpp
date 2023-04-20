//
// Created by mikibak on 4/20/23.
//

#include "Matrix.h"
#include <iostream>

Matrix::Matrix(int size_x, int size_y) {
    SIZE_X = size_x;
    SIZE_Y = size_y;

    columns = new int* [SIZE_X];
    for (int y = 0; y < SIZE_X; y++) {
        columns[y] = new int[SIZE_Y];
        for (int x = 0; x < SIZE_X; x++) {
            columns[y][x] = 0;
        }
    }
}

void Matrix::Print() {
    for(int y = 0; y < SIZE_Y; y++) {
        for(int x = 0; x < SIZE_X; x++) {
            std::cout << columns[y][x];
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

Matrix operator+(Matrix lhs, const Matrix& rhs) // otherwise, both parameters may be const references
{
    lhs += rhs; // reuse compound assignment
    return lhs; // return the result by value (uses move constructor)
}

Matrix& Matrix::operator+=(const Matrix& rhs){
    try {
        if(rhs.SIZE_X != this->SIZE_X && rhs.SIZE_Y != this->SIZE_Y) {
            throw std::runtime_error(std::string("Failed to add operands: incorrect width and height"));
        }

        if(rhs.SIZE_X != this->SIZE_X) {
            throw std::runtime_error(std::string("Failed to add operands: incorrect width, correct height"));
        }

        if(rhs.SIZE_Y != this->SIZE_Y) {
            throw std::runtime_error(std::string("Failed to add operands: correct width, incorrect height"));
        }
    } catch(std::runtime_error err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "Operands of size [" << this->SIZE_X << ", " << this->SIZE_Y << "] and " << rhs.SIZE_X << ", " << rhs.SIZE_Y << std::endl;
        return *this;
    }

    for(int y = 0; y < SIZE_Y; y++) {
        for(int x = 0; x < SIZE_X; x++) {
            this->columns[y][x] += rhs.columns[y][x];
        }
    }
    return *this;
}


Matrix::~Matrix() {
    for (int i = 0; i < SIZE_X; i++) {
        delete columns[i];
    }
    delete columns;
}