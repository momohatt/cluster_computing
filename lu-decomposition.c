#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define SIZE 4

void printArray(int a[SIZE][SIZE])
{
    int i, j;
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            printf("%d ", a[i][j]);
        }
        putchar('\n');
    }
    printf("---------\n");
}

int main (int argc, char *argv[])
{
    int A[SIZE][SIZE];
    int L[SIZE][SIZE];
    int U[SIZE][SIZE];
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            scanf("%d", &A[i][j]);
        }
    }
    printArray(A);
    for (int p = 0; p < SIZE; p++) {
        L[p][p] = A[p][p];
        for (int i = p + 1; i < SIZE; i++) {
            L[i][p] = A[i][p];
            L[p][i] = 0;
        }
        U[p][p] = 1;
        for (int i = p + 1; i < SIZE; i++) {
            U[p][i] = A[p][i] / A[p][p];
            U[i][p] = 0;
        }
        for (int i = p + 1; i < SIZE; i++) {
            for (int j = p + 1; j < SIZE; j++) {
                A[i][j] -= U[p][j] * L[i][p];
            }
        }
    }

    printArray(L);
    printArray(U);

    return 0;
}

