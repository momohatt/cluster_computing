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

double vecdot(double a[], double b[])
{
    int i;
    double sum;
    for (i = 0; i < N; i++)
        sum += a[i] * b[i];
    return sum;
}

void vecadd(double ans[], double a[], double b[])
{
    int i;
    if (ans == a) {
        for (i = 0; i < N; i++)
            a[i] += b[i];
        ans = a;
    } else if (ans == b) {
        for (i = 0; i < N; i++)
            b[i] += a[i];
        ans = b;
    } else {
        for (i = 0; i < N; i++)
            ans[i] = a[i] + b[i];
    }
}

void vecsub(double ans[], double a[], double b[])
{
    int i;
    if (ans == a) {
        for (i = 0; i < N; i++)
            a[i] -= b[i];
        ans = a;
    } else {
        for (i = 0; i < N; i++)
            ans[i] = a[i] - b[i];
    }
}

void vecscalar(double ans[], double s, double a[])
{
    int i;
    for (i = 0; i < N; i++)
        ans[i] = a[i] * s;
}

void matvec(double ans[], double mat[N][N], double vec[])
{
    int i, j;
    if (ans == vec) {
        double tmp[N];
        for (i = 0; i < N; i++) {
            tmp[i] = 0;
            for (j = 0; j < N; j++) {
                tmp[i] += mat[i][j] * vec[j];
            }
        }
        for (i = 0; i < N; i++) {
            ans[i] = tmp[i];
        }
    } else {
        for (i = 0; i < N; i++) {
            ans[i] = 0;
            for (j = 0; j < N; j++) {
                ans[i] += mat[i][j] * vec[j];
            }
        }
    }
}

void matadd(double ans[N][N], double A[N][N], double B[N][N])
{
    int i, j, k;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            ans[i][j] = A[i][j] + B[i][j];
}

void matscalar(double ans[N][N], double s, double A[N][N])
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            ans[i][j] = s * A[i][j];
}

void mat_diag_add(double ans[N][N], double A[N][N], double d[N])
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            ans[i][j] = A[i][j] + ((i == j)? d[i] : 0);
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

void printMat(const double array[N][N])
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

void printVec(const double vec[N], int id, char* name)
{
    int i;
    printf("%s ", name);
    printf("%2d: ", id);
    for (i = 0; i < N; i++) {
        printf("%2.4f ", vec[i]);
    }
    printf("\n");
}

