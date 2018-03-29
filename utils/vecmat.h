#ifndef _VECMAT_H_
#define _VECMAT_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define N 10

// to[] = from[];
void veccp(double to[], double from[]);

// to[] = val[];
void vecfill(double to[], double val);

void matfill(double to[N][N], double val);

bool vecIsZero(double vec[]);

// a[] * b[];
double vecdot(const double a[], const double b[]);

// ans[] = a[] + b[];
double* vecadd(double a[], double b[]);

// ans[] = a[] - b[];
double* vecsub(double a[], double b[]);

// ans = s * a[];
double* vecscalar(double s, double a[]);

double* matvec(double mat[N][N], double vec[]);

void matmat(double ans[N][N], double A[N][N], double B[N][N]);

// 転置
void transpose(double A[N][N]);

void printMat(double array[N][N]);
void printVec(double vec[], int id);


#endif // _VECMAT_H_
