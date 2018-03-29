#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define N 10

void printMat(double A[N][N + 1])
{
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N + 1; j++) {
            printf("%.1f ", A[i][j]);
        }
        putchar('\n');
    }
}

void gaussian_elimination(double A[N][N + 1], double x[])
{
    int i, j, k;

    // Forward elimination
    // - 対角要素をのぞいた左下要素を全て0にする(上三角行列にする)
    for (k = 0; k < N - 1; k++) {
        double akk = A[k][k];
        for (i = k + 1; i < N; i++) {
            double aik = A[i][k];
            for (j = k; j <= N; j++) {
                A[i][j] -= aik * (A[k][j] / akk);
            }
        }
    }
    printMat(A);

    // Back substitution
    x[N - 1] = A[N - 1][N] / A[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {
        double ax = 0.0;
        for (int j = i + 1; j < N; j++) {
            ax += A[i][j] * x[j];
        }
        x[i] = (A[i][N] - ax) / A[i][i];
    }
}

int main (int argc, char *argv[])
{
    double A[N][N + 1];
    double x[N];
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0 ; j < N; j++) {
            scanf("%lf", &A[i][j]);
        }
    }
    for (i = 0; i < N; i++) {
        scanf("%lf", &A[i][N]);
    }

    gaussian_elimination(A, x);

    for (i = 0; i < N; i++)
        printf("x[%d] = %.4f\n", i, x[i]);
    return 0;
}

