#ifndef MAT_DECOMP_H
#define MAT_DECOMP_H

#include <stdio.h>
#include <stdlib.h>

#define N 10

// L * x = bとなるxを求める
void forward_substitution(const double L[N][N], const double b[], double x[]);

// U * x = bとなるxを求める
void backward_substitution(const double U[N][N], const double b[], double x[]);

// L * U * u = rとなるuを求める　
void icres(const double L[N][N], const double U[N][N], const double r[N], double u[N]);

void gaussian_elimination(double A[N][N], double b[], double x[]);

// LU decomposition
void lu(double A[N][N], double L[N][N], double U[N][N]);

// Incomplete Cholesky decomposition

#endif // MAT_DECOMP_H
