//#include <stdio.h>
//#include <math.h>
//#include <time.h>
//#pragma warning(disable:4996)
//
////int a[STR][COLS] =
////{ {1, 4, 5},
//// {2, 6, 8},
//// {1, 0, 9},
//// {4, 2, 8} };
//
//void PrintMatrix(double* mass, int size) {
//	for (int str = 0; str < size; str++) {
//		for (int col = 0; col < size; col++) {
//			printf("%.16f ", mass[str * size + col]);
//		}
//		printf("\n");
//	}
//}
//
//void PrintVector(double* vect, int size) {
//	printf("\n");
//	for (int i = 0; i < size; i++) {
//		printf("%.15f ", vect[i]);
//	}
//
//}
//
//void Copy(double* from, double* to, int size) {
//	for (int i = 0; i < size; i++) {
//		for (int j = 0; j < size; j++) {
//			to[size * i + j] = from[size * i + j];
//		}
//	}
//}
//
//void CopyVector(double* newVec, double* b, int size) {
//	for (int i = 0; i < size; i++) {
//		newVec[i] = b[i];
//	}
//}
//
////void LU(double* L, double* U, double* A, int size) {
////	Copy(A, U, size);
////	for (int i = 0; i < size; i++)
////		for (int j = i; j < size; j++)
////			L[j*size+i] = U[j * size + i] / U[i * size + i];
////
////	for (int k = 1; k < size; k++)
////	{
////		for (int i = k - 1; i < size; i++)
////			for (int j = i; j < size; j++)
////				L[j*size+i] = U[j * size + i] / U[i * size + i];
////
////		for (int i = k; i < size; i++)
////			for (int j = k - 1; j < size; j++)
////				U[i*size+j] = U[i * size + j] - L[ i*size + (k - 1)] * U[(k - 1)*size+j] ;
////	}
////}
//
//void LU(double* A, double* L, double* U, int n)
//{
//	Copy(A, U, n);
//
//	for (int i = 0; i < n; i++)
//		for (int j = i; j < n; j++)
//			L[j * n + i] = U[j * n + i] / U[i * n + i];
//
//	for (int k = 1; k < n; k++)
//	{
//		for (int i = k - 1; i < n; i++)
//			for (int j = i; j < n; j++)
//				L[j * n + i] = U[j * n + i] / U[i * n + i];
//
//		for (int i = k; i < n; i++)
//			for (int j = k - 1; j < n; j++)
//				U[i * n + j] = U[i * n + j] - L[i * n + (k - 1)] * U[(k - 1) * n + j];
//	}
//
//}
//
//void ReadVector(char* filename, double* b, int size) {
//	FILE* fp = fopen(filename, "r");
//	if (fp == NULL) {
//		printf("\nFile hasn't been opened\n");
//		return;
//	}
//	for (int i = 0; i < size; i++) {
//		double elem;
//		fscanf(fp, "%lf ", &elem);
//		b[i] = elem;
//	}
//}
//
//void VectorInFile(FILE* fp, double* vector, int size) {
//	for (int i = 0; i < size; i++) {
//		fprintf(fp, "%.15f ", vector[i]);
//	}
//}
//
//void SolveEq(double* L, double* U, double* A, double* b, double* x, const int size) {
//	double* y = (double*)malloc(sizeof(double) * size);
//	double sum = 0;
//	y[0] = b[0];
//	for (int i = 1; i < size; i++) { // идем по строкам
//		sum = 0;
//		for (int k = 0; k <= i - 1; k++) {
//			sum = sum + L[i * size + k] * y[k];
//		}
//		y[i] = b[i] - sum;     // коэффициент при y[i]=1
//	}
//	int m = size - 1;
//	x[m] = y[m] / U[size * m + m];
//
//	for (int i = size - 2; i >= 0; i--) {
//		sum = 0;
//		for (int k = size - 1; k >= i + 1; k--) {
//			sum = sum + U[size * i + k] * x[k];
//		}
//		x[i] = (y[i] - sum) / U[size * i + i];
//	}
//	free(y);
//}
//
//void ReadMatrix(FILE* fp, double* A, int size) {
//	for (int i = 0; i < size; i++) {
//		for (int j = 0; j < size; j++) {
//			double elem = 0;
//			fscanf(fp, "%lf ", &elem);
//			A[i * size + j] = elem;
//		}
//	}
//}
//
//void FillZero(double* mass, int size) {
//	for (int i = 0; i < size; i++) {
//		for (int j = 0; j < size; j++) {
//			mass[i * size + j] = 0;
//		}
//	}
//}
//
//void Check(double* A, double* b, double* x, int size) {  //только первая строка
//	double sum = 0;
//	for (int i = 0; i < size; i++) {
//		sum = sum + A[i] * x[i];
//	}
//	printf("\n\n%lf %lf", b[0], sum);
//}
//
//void ChangeVector(double* b, double delta, int size) {
//	for (int i = 0; i < size; i++) {
//		b[i] = b[i] + delta;
//	}
//}
//
//int main() {
//	int size = 10;
//	int n = 9;
//	char* matrixFile = "C:/Users/z.kate/Desktop/3 сем/численные методы/matrix.csv";
//	char* vectorFile = "C:/Users/z.kate/Desktop/3 сем/численные методы/vector.csv";
//	char* xFile = "C:/Users/z.kate/Desktop/3 сем/численные методы/res.csv";
//	/*char* goodFile = "C:/Users/z.kate/Desktop/3 сем/chmData/лаба2/good.csv";
//	char* badFile = "C:/Users/z.kate/Desktop/3 сем/chmData/лаба2/bad.csv";
//	char* myMat = "C:/Users/z.kate/Desktop/3 сем/chmData/лаба2/myMat.csv";
//	char* newvect = "C:/Users/z.kate/Desktop/3 сем/chmData/2graph/newb.csv";
//	char* newresgood = "C:/Users/z.kate/Desktop/3 сем/chmData/2graph/newres.csv";
//	char* newresbad = "C:/Users/z.kate/Desktop/3 сем/chmData/2graph/newres2.csv";*/
//	FILE* fp1 = fopen(matrixFile, "r");
//	FILE* fp2 = fopen(vectorFile, "r");
//	FILE* fp3 = fopen(xFile, "w");
//	/*FILE* fpgood = fopen(goodFile, "w");
//	FILE* fpbad = fopen(badFile, "w");
//	FILE* fpnewb = fopen(newvect, "r");
//	FILE* fpnewresgood = fopen(newresgood, "w");
//	FILE* fpnewresbad = fopen(newresbad, "w");*/
//
//	double** massMat = (double**)malloc(sizeof(double*) * size);
//	double** massVec = (double**)malloc(sizeof(double*) * size);
//	double** massL = (double**)malloc(sizeof(double*) * size);
//	double** massU = (double**)malloc(sizeof(double*) * size);
//	double** massX = (double**)malloc(sizeof(double*) * size);
//	double** newb = (double**)malloc(sizeof(double*) * size);
//	double** newXgood = (double**)malloc(sizeof(double*) * n);
//	double** newXbad = (double**)malloc(sizeof(double*) * n);
//
//	
//
//	
//
//	//заполняем массив матриц
//	double* A;
//	for (int k = 0; k < 10; k++) {
//		A = (double*)malloc(sizeof(double) * size * size);
//		for (int i = 0; i < size; i++) {
//			for (int j = 0; j < size; j++) {
//				double elem = 0;
//				fscanf(fp1, "%lf ", &elem);
//				A[i * size + j] = elem;
//			}
//		}
//
//		massMat[k] = A;
//	}
//	double* B = massMat[0];
//	double* C = massMat[9];
//	
//
//	//заполняем массив векторов newb
//	double* b;
//	for (int k = 0; k < 10; k++) {
//		b = (double*)malloc(sizeof(double) * size);
//		for (int i = 0; i < size; i++) {
//			double elem;
//			fscanf(fp2, "%lf ", &elem);
//			b[i] = elem;
//		}
//		massVec[k] = b;
//	}
//	
//	
//	PrintVector(massVec[0], 10);
//	//заполним массив L U
//	double* L;
//	double* U;
//	for (int k = 0; k < 10; k++) {
//		L = (double*)malloc(sizeof(double) * size * size);
//		U = (double*)malloc(sizeof(double) * size * size);
//		FillZero(L, size);
//		FillZero(U, size);
//		LU(massMat[k], L, U, size);
//		massL[k] = L;
//		massU[k] = U;
//	}
//
//	//решаем слау и заполняем массив Х
//	double* x;
//	for (int k = 0; k < size; k++) {
//		x = (double*)malloc(sizeof(double) * size);
//		SolveEq(massL[k], massU[k], massMat[k], massVec[k], x, size);
//		massX[k] = x;
//	}
//
//	//запичываем результат х в файл
//	for (int k = 0; k < size; k++) {
//		for (int i = 0; i < size; i++) {
//			double* A = massX[k];
//			double elem = A[i];
//			fprintf(fp3, "%.16f ", elem);
//		}
//		fprintf(fp3, "\n");
//	}
//
//	////создаем два файла для 2 матриц +done
//	////берем 1 матрицу, в цикле вносим возмущение в вектор b, считаем х и записываем в файл
//	//
//	//// запись в файл для хорошй матрицы
//	//for (int k = 0; k < size; k++) {
//	//	//double* newVec = (double*)malloc(sizeof(double) * size);
//	//	double* newX = (double*)malloc(sizeof(double) * size);
//	//	double* Ltemp = (double*)malloc(sizeof(double) * size * size);
//	//	double* Utemp = (double*)malloc(sizeof(double) * size * size);
//	//	FillZero(Ltemp, size);
//	//	//CopyVector(newVec, b1, size);
//	//	LU(B, Ltemp, Utemp, size);
//	//	//ChangeVector(newVec, Deltab[k], size);
//	//	SolveEq(Ltemp, Utemp, B, newb[k], newX, size);
//	//	VectorInFile(fpnewresgood, newX, size);
//	//	fprintf(fpnewresgood, "\n");
//	//}
//
//	////// запись в файл для плохой матрицы
//	//for (int k = 0; k < size; k++) {
//	//	//double* newVec = (double*)malloc(sizeof(double) * size);
//	//	double* newX = (double*)malloc(sizeof(double) * size);
//	//	double* Ltemp = (double*)malloc(sizeof(double) * size * size);
//	//	double* Utemp = (double*)malloc(sizeof(double) * size * size);
//	//	FillZero(Ltemp, size);
//	//	//CopyVector(newVec, b2, size);
//	//	LU(C, Ltemp, Utemp, size);
//	//	//ChangeVector(newVec, Deltab[k], size);
//	//	SolveEq(Ltemp, Utemp, C, newb[k], newX, size);
//	//	PrintVector(newX, size);
//	//	VectorInFile(fpnewresbad, newX, size);
//	//	fprintf(fpnewresbad, "\n");
//	//}
//	//PrintMatrix(B, 10);
//	//closing all files
//	fclose(fp1);
//	fclose(fp2);
//	fclose(fp3);
//	/*fclose(fpgood);
//	fclose(fpbad);*/
//}