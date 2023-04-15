class Matrix {
	int SIZE_X;
	int SIZE_Y;
	int** columns;

	Matrix() {
		SIZE_X = 9 * 6 * 8;
		SIZE_Y = 9 * 6 * 8;

		columns = new int* [SIZE_X];
		for (int i = 0; i < SIZE_X; i++) {
			columns[i] = new int[SIZE_Y];
		}

		InitializeBandMatrix;
	}

	void InitializeBandMatrix() {

	}

	~Matrix() {
		for (int i = 0; i < SIZE_X; i++) {
			delete[] columns[i];
		}
		delete[] columns;
	}
};