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

void printMat(double array[N][N])
{
    int i, j;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%f ", array[i][j]);
        }
        putchar('\n');
    }
}

void printVec(double vec[N], int id)
{
    int i;
    printf("%2d: ", id);
    for (i = 0; i < N; i++) {
        printf("%.4f ", vec[i]);
    }
    putchar('\n');
}
