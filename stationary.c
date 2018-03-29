#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "utils/vecmat.h"

#define N 10
#define EPS  (1.0e-8)

bool hasConverged(double x[], double A[N][N], double b[])
{
    double ax[N], r[N]; // r = b - Ax
    matvec(ax, A, x);
    vecsub(r, b, ax);
    return (vecdot(r, r) < EPS)? true : false;
}

void gauss_seidel(double x[], double A[N][N], double b[])
{
    int i, j, k = 0;

    while (!hasConverged(x, A, b)) {
        k++;
        for (i = 0; i < N; i++) {
            x[i] = b[i];
            for (j = 0; j < N; j++) {
                if (i == j) continue;
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }
        printVec(x, k, "x");
    }
}

void sor(double x[], double A[N][N], double b[])
{
    double y[N] = {};
    double omega;
    int i, j, k = 0;

    printf("omega = ");
    scanf("%lf", &omega);

    while (!hasConverged(x, A, b)) {
        k++;
        for (i = 0; i < N; i++) {
            y[i] = b[i];
            for (j = 0; j < N; j++) {
                if (i == j) continue;
                y[i] -= A[i][j] * x[j];
            }
            y[i] /= A[i][i];
            x[i] += omega * (y[i] - x[i]);
        }
        printVec(x, k, "x");
    }
}

int main (int argc, char *argv[])
{
    double A[N][N];
    double b[N];
    double x[N] = {}; // 解
    vecfill(x, 0);
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            scanf("%lf", &A[i][j]);
    for (i = 0; i < N; i++)
        scanf("%lf", &b[i]);

    if (argc < 2) {
        printf("usage: ./[bin] [ sor | gs ]\n");
        return 0;
    }
    if (strcmp(argv[1], "sor")) {
        sor(x, A, b);
        printVec(x, -1, "x");
    } else if (strcmp(argv[1], "gs")) {
        gauss_seidel(x, A, b);
        printVec(x, -1, "x");
    }
    return 0;
}
