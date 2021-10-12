#include<stdio.h>



void PrintMatrix(int** mass, int size) {
	int col = 0;
	int str = 0;
	for (str = 0; str < size; str++) {
		for (col = 0; col < size; col++) {
			printf("%d ", *(mass+str*size+col));
		}
		printf("\n");
	}
}

int main() {
	int mass[3][3] = { {1, 2, 3}, {4, 5, 6}, {7, 8, 9} };
	int mass[5][6] = { {1, 2, 3} };
	int size = 3;
	PrintMatrix(&mass, size);

}