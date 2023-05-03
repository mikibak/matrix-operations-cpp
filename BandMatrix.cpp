#include "BandMatrix.h"

BandMatrix::BandMatrix(int sizeX, int sizeY, int a1, int a2, int a3) {
    this->ROWS = sizeY;
    this->COLS = sizeX;

    columns = new double* [ROWS];
    for(int y = 0; y < ROWS; y++) {
        columns[y] = new double[COLS];
        for(int x = 0; x < COLS; x++) {
            if(x == y) {
                //main diagonal
                columns[y][x] = a1;
            }
            if(x == y + 1 || y == x + 1) {
                //second diagonal
                columns[y][x] = a2;
            }
            if(x == y + 2 || y == x + 2) {
                //main diagonal
                columns[y][x] = a3;
            }
        }
    }
}
