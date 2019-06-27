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
    return vecnorm(r) < EPS;
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

void sor(double x[], double A[N][N], double b[], double omega)
{
    double y[N] = {};
    int i, j, k = 0;

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

void read_input(double A[N][N], double b[])
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            scanf("%lf", &A[i][j]);
    for (i = 0; i < N; i++)
        scanf("%lf", &b[i]);
}

int main (int argc, char *argv[])
{
    double A[N][N];
    double b[N];
    double x[N] = {}; // 解
    vecfill(x, 0);
    int i, j;

    if (argc < 2) {
        printf("usage: ./[bin] [ sor | gs ] [(omega)]\n");
        return 1;
    }

    if (strcmp(argv[1], "sor") == 0) {
        if (argc < 3) {
            printf("invalid input. usage: ./[bin] sor [omega]\n");
            return 0;
        }
        double omega = atof(argv[2]);
        if (omega <= 0 || omega >= 2) {
            printf("omega should be within range of (0, 2)\n");
            return 0;
        }

        printf("SOR method (omega = %f)\n", omega);
        read_input(A, b);
        sor(x, A, b, omega);
    } else if (strcmp(argv[1], "gs") == 0) {
        printf("Gauss Seidel method\n");
        read_input(A, b);
        gauss_seidel(x, A, b);
    } else {
        printf("invalid input. usage: ./[bin] [ sor | gs ]\n");
        return 1;
    }

    for (i = 0; i < N; i++) {
        printf("x[%d] = %f\n", i, x[i]);
    }
    return 0;
}
