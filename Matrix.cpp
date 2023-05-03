//
// Created by mikibak on 4/20/23.
//

#include "Matrix.h"
#include <iostream>
#include <cmath>

Matrix::Matrix(int cols, int rows) {
    COLS = cols;
    ROWS = rows;

    columns = new double* [ROWS];
    for (int row = 0; row < ROWS; row++) {
        columns[row] = new double[COLS];
        for (int col = 0; col < COLS; col++) {
            columns[row][col] = 0;
        }
    }
}

Matrix::Matrix(int cols, int rows, int value) {
    COLS = cols;
    ROWS = rows;

    columns = new double* [ROWS];
    for (int row = 0; row < ROWS; row++) {
        columns[row] = new double[COLS];
        for (int col = 0; col < COLS; col++) {
            columns[row][col] = value;
        }
    }
}

Matrix::Matrix() {
    COLS = 0;
    ROWS = 0;
}

void Matrix::Print() const{
    for(int y = 0; y < ROWS; y++) {
        for(int x = 0; x < COLS; x++) {
            std::cout << columns[y][x] << "   ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

Matrix operator+(const Matrix& lhs, const Matrix& rhs) // otherwise, both parameters may be const references
{
    Matrix* matrix = new Matrix(lhs.COLS, lhs.ROWS);

    matrix->COLS = lhs.COLS;
    matrix->ROWS = lhs.ROWS;

    for (int y = 0; y < matrix->ROWS; y++) {
        for (int x = 0; x < matrix->COLS; x++) {
            matrix->columns[y][x] = lhs.columns[y][x] + rhs.columns[y][x];
        }
    }
    return *matrix; // return the result by value (uses move constructor)
}

Matrix operator-(const Matrix& lhs, const Matrix& rhs) // otherwise, both parameters may be const references
{
    Matrix* matrix = new Matrix(lhs.COLS, lhs.ROWS);

    matrix->COLS = lhs.COLS;
    matrix->ROWS = lhs.ROWS;

    for (int y = 0; y < matrix->ROWS; y++) {
        for (int x = 0; x < matrix->COLS; x++) {
            matrix->columns[y][x] = lhs.columns[y][x] - rhs.columns[y][x];
        }
    }
    return *matrix; // return the result by value (uses move constructor)
}

Matrix& Matrix::operator+=(const Matrix& rhs){
    try {
        this->CheckIfSizeEqual(rhs);

    } catch(std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "Addition error with operands of size [" << this->COLS << ", " << this->ROWS << "] and [" << rhs.COLS << ", " << rhs.ROWS << "]" << std::endl;
        return *this;
    }

    for(int y = 0; y < ROWS; y++) {
        for(int x = 0; x < COLS; x++) {
            this->columns[y][x] += rhs.columns[y][x];
        }
    }
    return *this;
}

Matrix::~Matrix() {
    for (int i = 0; i < ROWS; i++) {
        delete columns[i];
    }
    delete columns;
}

double Matrix::GetElement(int col, int row) {
    try {
        this->CheckIfValidIndex(col, row);
    } catch(std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "Get element error with operands of size [" << this->COLS << ", " << this->ROWS << "] and indices [" << col << ", " << row << "]" << std::endl;
        return 0;
    }
    return this->columns[row][col];
}

void Matrix::SetElement(int col, int row, double value) {
    try {
        this->CheckIfValidIndex(col, row);
    } catch(std::runtime_error& err) {
        std::cerr << err.what() << std::endl;
        std::cerr << "Set element error with operands of size [" << this->COLS << ", " << this->ROWS << "] and indices [" << col << ", " << row << "]" << std::endl;
    }
    this->columns[row][col] = value;
}

void Matrix::CheckIfSizeEqual(const Matrix &other) {
    if(other.COLS != this->COLS && other.ROWS != this->ROWS) {
        throw std::runtime_error(std::string("Wrong operands: incorrect width and height"));
    }

    if(other.COLS != this->COLS) {
        throw std::runtime_error(std::string("Wrong operands: incorrect width, correct height"));
    }

    if(other.ROWS != this->ROWS) {
        throw std::runtime_error(std::string("Wrong operands: correct width, incorrect height"));
    }
}

void Matrix::CheckIfValidIndex(int x, int y) {
    if(x > this->COLS && y > this->ROWS) {
        throw std::runtime_error(std::string("Wrong indices: too large width and height"));
    }

    if(x > this->COLS) {
        throw std::runtime_error(std::string("Wrong index: too large width, correct height"));
    }

    if(y > this->ROWS) {
        throw std::runtime_error(std::string("Wrong index: correct width, too large height"));
    }

    if(x < 0 || y < 0) {
        throw std::runtime_error(std::string("Wrong indices: one or two indices below zero"));
    }
}

Matrix* Matrix::LowerTriangle() {
    Matrix* triangle = new Matrix(COLS, ROWS);

    for (int y = 0; y < ROWS; y++) {
        for (int x = 0; x < COLS; x++) {
            if(y > x) {
                triangle->SetElement(x, y, this->columns[y][x]);
            }
        }
    }

    return triangle;
}

Matrix* Matrix::UpperTriangle() {
    Matrix* triangle = new Matrix(COLS, ROWS);

    for (int y = 0; y < ROWS; y++) {
        for (int x = 0; x < COLS; x++) {
            if(y < x) {
                triangle->SetElement(x, y, this->columns[y][x]);
            }
        }
    }

    return triangle;
}

Matrix* Matrix::Diagonal() {
    Matrix* diagonal = new Matrix(COLS, ROWS);

    for (int y = 0; y < ROWS; y++) {
        for (int x = 0; x < COLS; x++) {
            if(y == x) {
                diagonal->SetElement(x, y, this->columns[y][x]);
            }
        }
    }

    return diagonal;
}

Matrix operator*(const Matrix& lhs, const Matrix& rhs) // otherwise, both parameters may be const references
{
    Matrix* matrix = new Matrix(rhs.COLS, lhs.ROWS);

    if(lhs.COLS != rhs.ROWS) {
        std::cerr << "Multiplication error with operands of size [" << lhs.ROWS << ", " << lhs.COLS << "] and [" << rhs.ROWS << ", " << rhs.COLS << "]" << std::endl;
        std::cerr << "Size Y of first operand must be the same as size X of the second one" << std::endl;
        return *matrix;
    }

    int sum;

    for(int i = 0; i < lhs.ROWS; i++) {
        for(int j = 0; j < rhs.COLS; j++) {
            sum = 0;
            for(int k = 0; k < lhs.COLS; k++) {
                sum += lhs.columns[i][k] * rhs.columns[k][j];
            }
            matrix->columns[i][j] = sum;
        }
    }

    return *matrix; // return the result by value (uses move constructor)
}

Matrix& Matrix::operator*=(const Matrix& rhs) {

    Matrix new_matrix = Matrix(this->COLS, rhs.ROWS);

    if(this->ROWS != rhs.COLS) {
        std::cerr << "Multiplication error with operands of size [" << this->COLS << ", " << this->ROWS << "] and [" << rhs.COLS << ", " << rhs.ROWS << "]" << std::endl;
        std::cerr << "Size Y of first operand must be the same as size X of the second one" << std::endl;
        return new_matrix;
    }

    int sum;

    for(int i = 0; i < COLS; i++) {
        for(int j = 0; j < rhs.ROWS; j++) {
            sum = 0;
            for(int k = 0; k < ROWS; k++) {
                sum += columns[i][k] * rhs.columns[k][j];
            }
            new_matrix.columns[i][j] = sum;
        }
    }

    this->COLS = new_matrix.COLS;
    this->ROWS = new_matrix.ROWS;
    for (int y = 0; y < ROWS; y++) {
        for (int x = 0; x < COLS; x++) {
            columns[y][x] = new_matrix.columns[y][x];
        }
    }

    return *this;
}

Matrix& Matrix::operator=(const Matrix& rhs) noexcept {

    this->COLS = rhs.COLS;
    this->ROWS = rhs.ROWS;
    for (int y = 0; y < ROWS; y++) {
        for (int x = 0; x < COLS; x++) {
            columns[y][x] = rhs.columns[y][x];
        }
    }

    return *this;
}

Matrix operator/(const Matrix& m1, const Matrix& m2) {
    if (m1.COLS != m2.COLS || m1.ROWS != m2.ROWS) {
        // handle error: matrices have different sizes
        throw std::invalid_argument("Matrices have different sizes.");
    }
    Matrix result(m1.COLS, m1.ROWS);
    for (int y = 0; y < m1.ROWS; y++) {
        for (int x = 0; x < m1.COLS; x++) {
            if (m2.columns[y][x] == 0) {
                // handle error: division by zero
                throw std::invalid_argument("Division by zero.");
            }
            result.columns[y][x] = m1.columns[y][x] / m2.columns[y][x];
        }
    }
    return result;
}

int Matrix::GetSizeY() const{
    return ROWS;
}

int Matrix::GetSizeX() const {
    return COLS;
}

double Matrix::Norm() const{
    double sum = 0.0;
    for (int y = 0; y < ROWS; y++) {
        for (int x = 0; x < COLS; x++) {
            sum += columns[y][x] * columns[y][x];
        }
    }
    return sqrt(sum);
}

Matrix Matrix::forwardSubstitution(const Matrix& b) const {
    if (COLS != ROWS) {
        throw std::invalid_argument("Matrix must be square");
    }
    Matrix x(1, ROWS, 1); // Initialize the solution vector to 1
    for (int i = 0; i < ROWS; i++) {
        double sum = 0;
        for (int j = 0; j < i; j++) {
            sum += this->columns[j][i] * x.columns[j][0]; // Calculate the sum of the previous solutions multiplied by the corresponding matrix coefficients
        }
        x.columns[i][0] = (b.columns[i][0] - sum) / this->columns[i][i]; // Calculate the ith solution using the updated sum and the ith diagonal coefficient
    }
    return x; // Return the solution vector
}
    /*Matrix* r = new Matrix(1, ROWS, 0);
    for (int row = 0; row < ROWS; row++) {
        int sum = 0;
        for (int column = 0; column < row; column++) {
            sum += columns[row][column] * r->columns[row][0];
        }
        int diagonal = columns[row][row];
        if (diagonal == 0) {
            throw std::invalid_argument("Matrix is not lower triangular");
        }
        r->columns[row][0] = (b.columns[row][0] - sum) / diagonal;
    }
    return *r;*/
//}