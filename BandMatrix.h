//
// Created by mikibak on 4/20/23.
//

#ifndef MATRIX_OPERATIONS_CPP_BANDMATRIX_H
#define MATRIX_OPERATIONS_CPP_BANDMATRIX_H


#include "Matrix.h"

class BandMatrix : public Matrix{
public:
    BandMatrix(int sizeX, int sizeY, int a1, int a2, int a3);
};


#endif //MATRIX_OPERATIONS_CPP_BANDMATRIX_H
