#include "vecmat.h"

void veccp(double to[], double from[])
{
    int i;
    for (i = 0; i < N; i++)
        to[i] = from[i];
}

void vecfill(double to[], double val)
{
    int i;
    for (i = 0; i < N; i++)
        to[i] = val;
}

void matfill(double to[N][N], double val)
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            to[i][j] = val;
}

double vecdot(const double a[], const double b[])
{
    int i;
    double sum;
    for (i = 0; i < N; i++)
        sum += a[i] * b[i];
    return sum;
}

double* vecadd(double a[], double b[])
{
    double *ans = malloc(sizeof(double) * N);
    int i;
//    if (ans == a) {
//        for (i = 0; i < N; i++)
//            a[i] += b[i];
//        ans = a;
//    } else if (ans == b) {
//        for (i = 0; i < N; i++)
//            b[i] += a[i];
//        ans = b;
//    } else {
        for (i = 0; i < N; i++)
            ans[i] = a[i] + b[i];
//    }
    return ans;
}

double* vecsub(double a[], double b[])
{
    double *ans = malloc(sizeof(double) * N);
    int i;
//    if (ans == a) {
//        for (i = 0; i < N; i++)
//            a[i] -= b[i];
//        ans = a;
//    } else {
        for (i = 0; i < N; i++)
            ans[i] = a[i] - b[i];
//    }
    return ans;
}

double* vecscalar(double s, double a[])
{
    double *ans = malloc(sizeof(double) * N);
    int i;
    for (i = 0; i < N; i++)
        ans[i] = a[i] * s;
    return ans;
}

double* matvec(double mat[N][N], double vec[])
{
    double *ans = malloc(sizeof(double) * N);
    int i, j;
//    if (ans == vec) {
//        double tmp[N];
//        for (i = 0; i < N; i++) {
//            tmp[i] = 0;
//            for (j = 0; j < N; j++) {
//                tmp[i] += mat[i][j] * vec[j];
//            }
//        }
//        for (i = 0; i < N; i++) {
//            ans[i] = tmp[i];
//        }
//    } else {
        for (i = 0; i < N; i++) {
            ans[i] = 0;
            for (j = 0; j < N; j++) {
                ans[i] += mat[i][j] * vec[j];
            }
        }
//    }
    return ans;
}

void matmat(double ans[N][N], double A[N][N], double B[N][N])
{
    int i, j, k;
    if (ans == A) {
        double tmp[N][N];
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                tmp[i][j] = 0.0;
                for (k = 0; k < N; k++) {
                    tmp[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++)
                ans[i][j] = tmp[i][j];
    } else {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                ans[i][j] = 0.0;
                for (k = 0; k < N; k++) {
                    ans[i][j] += A[i][k] * B[k][j];
                }
            }
        }
    }
}

void transpose(double A[N][N])
{
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = i + 1; j < N; j++) {
            double tmp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = tmp;
        }
    }
}

void printMat(double array[N][N])
{
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%.4f ", array[i][j]);
        }
        putchar('\n');
    }
    putchar('\n');
}

void printVec(double vec[N], int id)
{
    int i;
    printf("%2d: ", id);
    for (i = 0; i < N; i++) {
        printf("%.4f ", vec[i]);
    }
    printf("\n\n");
}

