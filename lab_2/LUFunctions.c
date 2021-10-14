#include <stdio.h>
#include <math.h>
#pragma warning(disable:4996)

//int a[STR][COLS] =
//{ {1, 4, 5},
// {2, 6, 8},
// {1, 0, 9},
// {4, 2, 8} };


void PrintMatrix(double* mass, int size) {
	for (int str = 0; str < size; str++) {
		for (int col = 0; col < size; col++) {
			printf("%.8f ", mass[str*size+col]);
		}
		printf("\n");
	}
}

void PrintVector(double* vect, int size) {
	printf("\n");
	for (int i = 0; i < size; i++) {
		printf("%.8f\n", vect[i]);
	}

}

void Copy(double* from, double* to, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			to[size * i + j] = from[size * i + j];
		}
	}
}


//void LU(double* L, double* U, double* A, int size) {
//	Copy(A, U, size);
//	for (int i = 0; i < size; i++)
//		for (int j = i; j < size; j++)
//			L[j*size+i] = U[j * size + i] / U[i * size + i];
//
//	for (int k = 1; k < size; k++)
//	{
//		for (int i = k - 1; i < size; i++)
//			for (int j = i; j < size; j++)
//				L[j*size+i] = U[j * size + i] / U[i * size + i];
//
//		for (int i = k; i < size; i++)
//			for (int j = k - 1; j < size; j++)
//				U[i*size+j] = U[i * size + j] - L[ i*size + (k - 1)] * U[(k - 1)*size+j] ;
//	}
//}

void LU(double* A, double* L, double* U, int n)
{
	Copy(A, U, n);

	for (int i = 0; i < n; i++)
		for (int j = i; j < n; j++)
			L[j*n+i] = U[j*n+i] / U[i*n+i];

	for (int k = 1; k < n; k++)
	{
		for (int i = k - 1; i < n; i++)
			for (int j = i; j < n; j++)
				L[j*n+i] = U[j*n+i] / U[i*n+i];

		for (int i = k; i < n; i++)
			for (int j = k - 1; j < n; j++)
				U[i*n+j] = U[i*n+j] - L[i*n+(k - 1)] * U[(k - 1)*n+j];
	}

}


void ReadVector(char* filename, double* b, int size) {
	FILE* fp = fopen(filename, "r");
	if (fp == NULL) {
		printf("\nFile hasn't been opened\n");
		return;
	}
	for (int i = 0; i < size; i++) {
		double elem;
		fscanf(fp, "%lf ", &elem);
		b[i] = elem;
	}
}


void VectorInFile(char* filename, double* vector, int size) {
	FILE* fp = fopen(filename, "a");
	if (fp == NULL) {
		printf("file wan't opened");
		return;
	}
	for (int i = 0; i < size; i++) {
		fprintf(fp, "%.8f ", vector[i]);
	}
	fclose(fp);
}



void SolveEq(double* L, double* U, double* A, double* b, double* x, const int size) {
	double* y = (double*)malloc(sizeof(double) * size);
	double sum = 0;
	y[0] = b[0];
	for (int i = 1; i < size; i++) { // идем по строкам
		sum = 0;
		for (int k = 0; k <= i - 1; k++) {
			sum = sum + L[i * size + k] * y[k];
		}
		y[i] = b[i] - sum;     // коэффициент при y[i]=1
	}
	int m = size - 1;
	x[m] = y[m] / U[size * m + m];

	for (int i = size - 2; i >= 0; i--) {
		sum = 0;
		for (int k = size - 1; k >= i + 1; k--) {
			sum = sum + U[size * i + k] * x[k];
		}
		x[i] = (y[i] - sum) / U[size * i + i];
	}
	free(y);
}

void ReadMatrix(FILE* fp, double* A, int size) {
	for (int i = 0; i < size; i++) {
		for (int j = 0; j < size; j++) {
			double elem = 0;
			fscanf(fp, "%lf ", &elem);
			A[i * size + j] = elem;
		}
	}
}

void Print(double mass[][10], int size) {
	for (int str = 0; str < size; str++) {
		for (int col = 0; col < size; col++) {
			printf("%.8f ", *mass[str*size+col]);
		}
		printf("\n");
	}
}

int main() {
	int size = 10;
	double** massMat = (double**)malloc(sizeof(double*) * size );	
	double* massVec = (double*)malloc(sizeof(double) * size);
	
	double b[10] ;
	double L[10][10] = {0};
	double U[10][10] = {0};
	
	double* x = (double*)malloc(sizeof(double) * size);
	char* matrixFile = "C:/Users/z.kate/Desktop/3 сем/chmData/лаба2/matrix.csv";
	char* vectorFile = "C:/Users/z.kate/Desktop/3 сем/chmData/лаба2/vector.csv";
	FILE* fp1 = fopen(matrixFile, "r");
	FILE* fp2 = fopen(vectorFile, "r");
	
	//заполняем массив матриц
	double* A;
	for (int k = 0; k < 10; k++) {
		A = (double*)malloc(sizeof(double) * size*size);
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				double elem = 0;
				fscanf(fp1, "%lf ", &elem);
				A[i*size+j] = elem;
			}
		}
		massMat[k] = A;
	}

	
	PrintMatrix(massMat[1], size);
	LU(massMat[0], L, U, size);
	PrintMatrix(L, size);
	printf("\n");
	PrintMatrix(U, size);
	//PrintMatrix(A, size);
	//ReadVector(vectorFile, b, size);
	////PrintVector(b, size);
	//LU(A, L, U, size);
	//SolveEq(L, U, A, b, x, size);
	//PrintVector(x, size);
}

