#ifndef MAT_DECOMP_H
#define MAT_DECOMP_H

#include <stdio.h>
#include <stdlib.h>

#define N 10

void gaussian_elimination(double A[N][N + 1], double x[]);

// LU decomposition
void lu(double A[N][N], double L[N][N], double U[N][N]);

// Incomplete Cholesky decomposition

#endif // MAT_DECOMP_H
