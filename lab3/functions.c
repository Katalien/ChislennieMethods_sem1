#include <stdio.h>

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
    }
    return res;
}

void PrintMatrix(int** mass, int size) {
    printf("\n");
    for (int str = 0; str < size; str++) {
        for (int col = 0; col < size; col++) {
            printf("%d ", mass[str][col]);
            // printf("%d ", mass[str] [ col]);
        }
        printf("\n");
    }
}



int main() {
    int m = 2;
    int size = 2;
    double** A = CreateMatrix(size);
    double* b = (double*)malloc(sizeof(double) * size);
    A[0][0] = 1.15;
    A[0][1] = 0.53;
    A[1][0] = 0.53;
    A[1][1] = 2.84;
    b[0] = 1.69;
    b[1] = 3.38;

    return 0;
}