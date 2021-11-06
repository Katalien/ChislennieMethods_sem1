#define _USE_MATH_DEFINES
#include <stdio.h>
#include<math.h>


int** CreateMatrix(int n) {
    int** matrix = (int**)malloc(sizeof(int*) * n);
    if (matrix == NULL) {
        free(matrix);
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        matrix[i] = (int*)malloc(sizeof(int) * n);
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

void FillZero(int** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[i][j] = 0;
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

void PrintVector(double* vec, int size) {
    printf("\n");
    for (int i = 0; i < size; i++) {
        printf("%.3f ", vec[i]);
    }
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

double* VecMinVec(double* a, double* b, int size) {
    double* res = (double*)malloc(sizeof(double) * size);
    if (res == NULL) {
        return 0;
    }
    for (int i = 0; i < size; i++) {
        res[i] = a[i] - b[i];
        //printf("\n%f ", res[i]);
    }
    return res;
}

void CountT(double* t, int m) {
    for (int i = 1; i <= m; i++) {
        t[i-1] = cos(M_PI * (double)((2 * i - 1)) / (double)(2 * m));
       // printf("\n %f ", t[i-1]);
    }
}

void CountLambda(double* lam, double* t, double d1, double d2, int m) {
    for (int i = 0; i < m; i++) {
        lam[i] = (d1 + d2) / 2 + (d2 - d1) * t[i] / 2;
       // printf("\n%f ", lam[i]);
    }
}

void CountAlpha(double* alpha, double* lam, int m) {
    for (int i = 0; i < m; i++) {
        alpha[i] = 1 / lam[i];
        //printf("\n%f ", alpha[i]);
    }
}

double* g(double* x, double alpha, double** A, double* b, int size) {
    double* tmp = MatMulVec(A, x, size);
    tmp = VecMinVec(tmp, b, size );
    PrintVector(tmp, size);
    for (int i = 0; i < size; i++) {
        tmp[i] = alpha * tmp[i];
        //printf("\n%f ", tmp[i]);
    }
    double* res = VecMinVec(x, tmp, size);
    PrintVector(res, size);
    return res;
}

double Norma(double* x, int size) {
    double res = 0, sum = 0;
    for (int i = 0; i < size; i++) {
        sum += x[i] * x[i];
    }
    res = sqrt(sum);
    return res;
}

double* Richardson(double** A, double* b, double*x0, int size, double d1, double d2, int m, double eps) {
    double* x = (double*)malloc(sizeof(double) * size);
    double* tmp = (double*)malloc(sizeof(double) * size);
    double* xnext = (double*)malloc(sizeof(double) * size); //?
    double* t = (double*)malloc(sizeof(double) * m);
    double* lam = (double*)malloc(sizeof(double) * m);
    double* alpha = (double*)malloc(sizeof(double) * m);
    CountT(t, m);
    CountLambda(lam, t, d1, d2, m);
    CountAlpha(alpha, lam, m);
    //x = g(x0, alpha[0], A, b, size);
   // double* dif = VecMinVec(x, x0, size);
    double* dif;
    tmp = x0;
   do {
       int count = 0;
        for (int i = 0; i < m; i++) {
            xnext = g(tmp, alpha[i], A, b, size);
            tmp = xnext;
        }
        dif = VecMinVec(xnext, x0, size);
        printf("\ndef = %f\n", Norma(dif, size));
        printf("\n%d : x = ", count++);
        PrintVector(xnext, size);
        x0 = xnext;
   } while (Norma(dif, size) > eps);
    return xnext;
}

int main() {
    int m = 2;
    int size = 2;
    double** A = CreateMatrix(size);
    double* b = (double*)malloc(sizeof(double) * size);
    double* x0 = (double*)malloc(sizeof(double) * size);
    double d1 = 0.99-0.01;
    double d2 = 2.99+0.01;
    A[0][0] = 1.15;
    A[0][1] = 0.53;
    A[1][0] = 0.53;
    A[1][1] = 2.84;
    b[0] = 1.69;
    b[1] = 3.38;
    x0[0] = 0;
    x0[1] = 0;
    double* res = Richardson(A, b,x0, size, d1, d2, m, 0.000001);
    PrintVector(res, size);
   

    return 0;
}