#include "BandMatrix.h"

BandMatrix::BandMatrix(int sizeX, int sizeY, int a1, int a2, int a3) : Matrix(sizeX, sizeY) {
    for(int y = 0; y < SIZE_Y; y++) {
        for(int x = 0; x < SIZE_X; x++) {
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
