//
// Created by mikibak on 4/20/23.
//

#include "Matrix.h"
#include <iostream>

Matrix::Matrix(int size_x, int size_y) {
    SIZE_X = size_x;
    SIZE_Y = size_y;

    columns = new int* [SIZE_Y];
    for (int y = 0; y < SIZE_Y; y++) {
        columns[y] = new int[SIZE_X];
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
        this->CheckIfSizeEqual(rhs);

    } catch(std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "Addition error with operands of size [" << this->SIZE_X << ", " << this->SIZE_Y << "] and [" << rhs.SIZE_X << ", " << rhs.SIZE_Y << "]"<< std::endl;
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

int Matrix::GetElement(int x, int y) {
    try {
        this->CheckIfValidIndex(x, y);
    } catch(std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "Get element error with operands of size [" << this->SIZE_X << ", " << this->SIZE_Y << "] and indices [" << x << ", " << y << "]"<< std::endl;
        return 0;
    }
    return this->columns[y][x];
}

void Matrix::SetElement(int x, int y, int value) {
    try {
        this->CheckIfValidIndex(x, y);
    } catch(std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "Set element error with operands of size [" << this->SIZE_X << ", " << this->SIZE_Y << "] and indices [" << x << ", " << y << "]"<< std::endl;
    }
    this->columns[y][x] = value;
}

void Matrix::CheckIfSizeEqual(const Matrix &other) {
    if(other.SIZE_X != this->SIZE_X && other.SIZE_Y != this->SIZE_Y) {
        throw std::runtime_error(std::string("Wrong operands: incorrect width and height"));
    }

    if(other.SIZE_X != this->SIZE_X) {
        throw std::runtime_error(std::string("Wrong operands: incorrect width, correct height"));
    }

    if(other.SIZE_Y != this->SIZE_Y) {
        throw std::runtime_error(std::string("Wrong operands: correct width, incorrect height"));
    }
}

void Matrix::CheckIfValidIndex(int x, int y) {
    if(x > this->SIZE_X && y > this->SIZE_Y) {
        throw std::runtime_error(std::string("Wrong indices: too large width and height"));
    }

    if(x > this->SIZE_X) {
        throw std::runtime_error(std::string("Wrong index: too large width, correct height"));
    }

    if(y > this->SIZE_Y) {
        throw std::runtime_error(std::string("Wrong index: correct width, too large height"));
    }

    if(x < 0 || y < 0) {
        throw std::runtime_error(std::string("Wrong indices: one or two indices below zero"));
    }
}
