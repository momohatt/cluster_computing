#include "mat_decomp.h"

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

void lu(double A[N][N], double L[N][N], double U[N][N])
{
    int i, j;
    int p = 0;
    while (p < N) {
        L[p][p] = A[p][p];
        for (i = p + 1; i < N; i++) {
            L[i][p] = A[i][p];
            L[p][i] = 0;
        }
        U[p][p] = 1;
        for (i = p + 1; i < N; i++) {
            U[p][i] = A[p][i] / A[p][p];
            U[i][p] = 0;
        }
        for (i = p + 1; i < N; i++) {
            for (j = p + 1; j < N; j++) {
                A[i][j] -= U[p][j] * L[i][p];
            }
        }
        p++;
    }
}
