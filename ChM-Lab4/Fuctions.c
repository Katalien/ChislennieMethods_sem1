#define _USE_MATH_DEFINES
#include <stdio.h>
#include<math.h>
#include<stdlib.h>
#pragma warning(disable:4996)

#define WORK 1 // file

double** CreateMatrix(double n) {
    double** matrix = (double**)malloc(sizeof(double*) * n);
    if (matrix == NULL) {
        free(matrix);
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        matrix[i] = (double*)malloc(sizeof(double) * n);
        if (matrix[i] == NULL) {
            for (int j = 0; j < i; j++) {
                free(matrix[j]);
            }
            free(matrix);
            return NULL;
        }
    }
    return matrix;
}

void FillZero(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0.0;
        }
    }
}

void PrintMatrix(double** mass, int size) {
    printf("\n");
    for (int str = 0; str < size; str++) {
        for (int col = 0; col < size; col++) {
            printf("%.3f ", mass[str][col]);
        }
        printf("\n");
    }
}

void ReadMatrix(FILE* fp, double** A, int size) {
    double elem = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            fscanf(fp, "%lf", &elem);
            A[i][j] = elem;
        }
    }
}

void ReadVector(FILE* fp, double* b, int size) {
    double elem = 0;
    for (int i = 0; i < size; i++) {
        fscanf(fp, "%lf", &elem);
        b[i] = elem;
    }
}

void PrintVector(double* vec, int size) {
    printf("\n");
    for (int i = 0; i < size; i++) {
        printf("%.10f ", vec[i]);
    }
}

void FillVectorZero(double* b, int size) {
    for (int i = 0; i < size; i++) {
        b[i] = 0;
    }
}

void FillVectorOnes(double* b, int size) {
    for (int i = 0; i < size; i++) {
        b[i] = 1;
    }
}


double FindMax(double* b, int size) {
    double max = b[0];
    for (int i = 1; i < size; i++) {
        if (b[i] >= max) {
            max = b[i];
        }
    }
    return max;
}

double FindMin(double* b, int size) {
    double min = b[0];
    for (int i = 1; i < size; i++) {
        if (b[i] <= min) {
            min = b[i];
        }
    }
    return min;
}

double* MatMulVec(double** A, double* b, int size) {
    double* res = (double*)malloc(sizeof(double) * size);
    if (res == NULL) {
        return 0;
    }
    for (int i = 0; i < size; i++) {
        double sum = 0;
        for (int j = 0; j < size; j++) {
            sum += A[i][j] * b[j];
        }
        res[i] = sum;
        //printf("\n%f ", res[i]);
    }
    return res;
}

double MatrixNorm(double** A, int size) {
    double sum = 0;
    double norm = 0;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            sum = sum+ powf(A[i][j], 2);
        }
    }
    norm = sqrt(sum);
    return norm;
}

double* VecDivNum(double* b, double a, int size) {
    double* res = (double*)malloc(sizeof(double) * size);
    if (res == NULL) {
        printf("error in malloc");
        return NULL;
    }
    for (int i = 0; i < size; i++) {
        res[i] = b[i] / a;
    }
    return res;
}

double* VecMulNum(double* b, double num, int size) {
    double* res = (double*)malloc(sizeof(double) * size);
    if (res == NULL) {
        return 0;
    }
    for (int i = 0; i < size; i++) {
        res[i] = b[i] * num;
    }
    return res;
}

void MatMinEMat(double** A, double** B, double koef, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            B[i][j] = A[i][j];
            if (i == j) {
                B[i][j] = A[i][j] - koef;
            }
        }
    }
}

double PowerMethod(double** A, double* x0, double** lam, double eps, int size) {
    double* x = (double*)malloc(sizeof(double) * size);
    double* x_next = (double*)malloc(sizeof(double) * size);
    double d1, d2;
    int i = 1;
    x = x0;
    x_next = MatMulVec(A, x, size);
    d2 = x_next[0] / x[0];
    do {
        i++;
        d1 = d2;
        x = x_next;
        x_next = MatMulVec(A, x, size);
        d2 = x_next[0] / x[0];
    } while (fabs(d2 - d1) > eps);
    *lam= VecDivNum(x_next, pow(d2, i), size);
    PrintVector(*lam, size);
    return d2;
}

double ShiftPowerMethod(double** A, double* x0, double** lam, double eps, int size) {
    double** B = CreateMatrix(size);
    double m = MatrixNorm(A, size);
    double d_tmp = 0, d = 0;
    MatMinEMat(A, B, m, size);
    FillVectorZero(*lam, size);
    d_tmp = PowerMethod(B, x0, lam, eps, size);
    d = d_tmp + m;
    return d;
}

void CheckMethod(double** A, double d, double* lam, double eps, int size) {
    double* a = MatMulVec(A, lam, size);
    double* b = VecMulNum(lam, d, size);
    int bool = 0;
    for (int i = 0; i < size; i++) {
        printf("\n%.20lf", (fabs(a[i] - b[i])));
        if (fabs(a[i] - b[i]) > eps) {
            bool = 1;
        }
    }
    if (bool == 0) {
        printf("\n\nsuccess");
    }
    else {
        printf("\n\nfail");
    }
}

int main() {
#ifdef WORK {
    FILE* fp1 = fopen("C:/Users/z.kate/Desktop/3 сем/численные методы/lab4/matrix.csv", "r");
    int size = 5;
    double** A = CreateMatrix(size);
    double* x0 = (double*)malloc(sizeof(double) * size);
    double* lam = (double*)malloc(sizeof(double) * size); ;
    double eps = 1e-20;
    double d = 0;
    ReadMatrix(fp1, A, size);
    FillVectorOnes(x0, size);
    d = ShiftPowerMethod(A, x0, &lam, eps, size);
    printf("\nd = %lf\n", d);
    CheckMethod(A, d, lam, eps, size);
    free(A);
    free(x0);
    free(lam);
    return 0;
#endif
#ifndef  WORK
    int size = 2;
    double eps = 0.1;
    double** A = CreateMatrix(size);
    double* b = (double*)malloc(sizeof(double) * size);
    double* x0 = (double*)malloc(sizeof(double) * size);
    double* lam = (double*)malloc(sizeof(double) * size);
    double d = 0;
    A[0][0] = 2.49;
    A[0][1] = 5.10;
    A[1][0] = 5.10;
    A[1][1] = 18.5;
    b[0] = 1.69;
    b[1] = 3.38;
    x0[0] = 1;
    x0[1] = 1;
    d = ShiftPowerMethod(A, x0, &lam, eps, size);
    CheckMethod(A, d, lam,eps, size);
#endif // ! WORK_FILE
    }