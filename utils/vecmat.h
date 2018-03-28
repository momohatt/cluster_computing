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
double vecdot(double a[], double b[]);

// ans[] = a[] + b[];
void vecadd(double ans[], double a[], double b[]);

// ans[] = a[] - b[];
void vecsub(double ans[], double a[], double b[]);

// ans = s * a[];
void vecscalar(double ans[], double s, double a[]);

void matvec(double ans[], double mat[N][N], double vec[]);

void matmat(double ans[N][N], double A[N][N], double B[N][N]);

// 転置
void transpose(double A[N][N]);

void printMat(double array[N][N]);
void printVec(double vec[], int id);


#endif // _VECMAT_H_
