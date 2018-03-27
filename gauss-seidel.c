#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "utils/vecmat.h"

#define N 10
#define MAX 100
#define EPS  (1.0e-8)

bool hasConverged(double x[], double A[N][N], double b[])
{
    double ax[N], r[N]; // r = b - Ax
    matvec(ax, A, x);
    vecsub(r, b, ax);
    return (vecdot(r, r) < EPS)? true : false;
}

int main (int argc, char *argv[])
{
    double A[N][N] =
        { {5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0, 0.0},
          {0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0, 0.0},
          {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0, 0.0},
          {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0, 2.0},
          {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 5.0} };
    double b[N] = {3.0, 1.0, 4.0, 0.0, 5.0, -1.0, 6.0, -2.0, 7.0, -15.0};

    double x[N] = {};
    int i, j, k = 0;

    // TODO: terminating condition
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
        printVec(x, k);
    }

    return 0;
}
