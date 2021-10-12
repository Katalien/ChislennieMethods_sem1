#include<stdio.h>
//_CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

//int a[STR][COLS] =
//{ {1, 4, 5},
// {2, 6, 8},
// {1, 0, 9},
// {4, 2, 8} };


void PrintMatrix(double* mass, int size) {
	for (int str = 0; str < size; str++) {
		for (int col = 0; col < size; col++) {
			printf("%.2f ", mass[str*size+col]);
		}
		printf("\n");
	}
}

void PrintVector(double* vect, int size) {
	printf("\n");
	for (int i = 0; i < size; i++) {
		printf("%.2f\n", vect[i]);
	}

}

void Copy(double* from, double* to, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			to[size * i + j] = from[size * i + j];
		}
	}
}


void LU(double* L, double* U, double* A, int size) {
	Copy(A, U, size);
	for (int i = 0; i < size; i++)
		for (int j = i; j < size; j++)
			L[j*size+i] = U[j * size + i] / U[i * size + i];

	for (int k = 1; k < size; k++)
	{
		for (int i = k - 1; i < size; i++)
			for (int j = i; j < size; j++)
				L[j*size+i] = U[j * size + i] / U[i * size + i];

		for (int i = k; i < size; i++)
			for (int j = k - 1; j < size; j++)
				U[i*size+j] = U[i * size + j] - L[ i*size + (k - 1)] * U[(k - 1)*size+j] ;
	}
}

void SolveEq(double* L, double* U, double* A, double* b, double* x, const int size) {
	double* y = (double*)malloc(sizeof(double)*size);
	int sum = 0;
	y[0] = b[0];

	//обратная подстановка: L диагональ 1 
	for (int i = 1; i < size; i++) { // идем по строкам
		sum = 0;
		for (int k = 0; k <= i-1 ; k++) {
				sum = sum + L[i * size + k] * y[k];      // считаем сумму по сроке: коэффициент из матрицы L * нйденный у
		}
		y[i] = b[i] - sum;     // коэффициент при y[i]=1
	}
	int m = size - 1;
	x[m] = y[m] / U[size * m + m];

	for (int i = size - 2; i >= 0; i--) {
		sum = 0;
		for (int k = size - 1; k >= i+1; k--) {
			sum = sum + U[size * i + k] * x[k];
		}
		x[i] = (y[i] - sum) / U[size * i + i];
	}
	free(y);
}

void ReadMatrix(char*filename, double* A, int size) {
	FILE* fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("\nFile hasn't been opened\n");
		return;
	}
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			int elem;
			fscanf(fp, "%d ", &elem );
			A[i * size + j] = elem;
		}
	}
}

void VectorInFile(char* filename, double* vector, int size) {
	FILE* fp = fopen(filename, "a");
	if (fp == NULL) {
		printf("file wan't opened");
		return;
	}
	for (int i = 0; i < size; i++) {
		fprintf(fp, "%.2f ", vector[i]);
	}
	fclose(fp);
}

int main() {
	/*double A[3][3] = { {1, 1,1}, {1, 2, 3}, {1, 3,4} };
	double L[3][3] = { 0 };
	double U[3][3] = { 0 };
	double b[3] = { 6, 14, 19 };
	const int size = 3;
	double* x = (double*)malloc(sizeof(double)*size);*/
	
	double A[3][3];
	double b[3] = {10, 10};
	double L[3][3];
	double U[3][3];
	int size = 3;
	double* x = (double*)malloc(sizeof(double) * size);
	
	//VectorInFile("vector.txt", x, size);
	char* matrixFile = "C:/Users/z.kate/Desktop/matrix.txt";
	ReadMatrix(matrixFile, A, size );
	PrintMatrix(A, size);
}

