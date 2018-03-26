#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#define VECSIZE 10
#define MAX 100

void veccp(int to[], int from[])
{
    int i;
    for (i = 0; i < VECSIZE; i++)
        to[i] = from[i];
}

void vecfill(int to[], int val)
{
    int i;
    for (i = 0; i < VECSIZE; i++)
        to[i] = val;
}

bool vecIsZero(int vec[])
{
    int i;
    for (i = 0; i < VECSIZE; i++)
        if (vec[i]) return false;
    return true;
}

int vecdot(int a[], int b[])
{
    int sum, i;
    for (i = 0; i < VECSIZE; i++)
        sum += a[i] * b[i];
    return sum;
}

void vecadd(int ans[], int a[], int b[])
{
    int i;
    for (i = 0; i < VECSIZE; i++)
        ans[i] = a[i] + b[i];
}

void vecscalar(int ans[], int s, int a[])
{
    int i;
    for (i = 0; i < VECSIZE; i++)
        ans[i] = a[i] * s;
}

void matvec(int ans[], int mat[VECSIZE][VECSIZE], int vec[])
{
    int i, j;
    for (i = 0; i < VECSIZE; i++)
        for (j = 0; j < VECSIZE; j++)
            ans[i] = vecdot(mat[i], vec);
}

int main (int argc, char *argv[])
{
    int b[VECSIZE];
    int X[VECSIZE];
    int k = 0; // index
    int x[MAX][VECSIZE];
    int r[MAX][VECSIZE];
    int s[MAX][VECSIZE];
    int A[VECSIZE][VECSIZE];
    vecfill(x[0], 0);
    veccp(r[0], b);

    // preconditioning
    int M[VECSIZE][VECSIZE];
    int i, j;
    for (i = 0; i < VECSIZE; i++) {
        for (j = 0; j < VECSIZE; j++) {
            M[i][j] = (i == j)? (double) (1 / A[i][j]) : 0;
        }
    }
    for (i = 0; i < VECSIZE; i++) {
        for (j = 0; j < VECSIZE; j++) {
            A[i][j] = M[i][i] * A[i][j];
        }
    }
    matvec(b, M, b);


    while (!vecIsZero(r[k])) {
        k++;
        if (k == 1) {
            veccp(s[1], r[0]);
        } else {
            int beta = vecdot(r[k - 1], r[k - 1]) / vecdot(r[k - 2], r[k - 2]);
            int tmp[VECSIZE];
            vecscalar(tmp, beta, s[k - 1]);
            vecadd(s[k], r[k - 1], tmp);
        }
        int tmp1[VECSIZE], tmp2[VECSIZE];

        // alpha = (r[k-1] * r[k-1]) / (s[k] * A * s[k])
        matvec(tmp1, A, s[k]);
        int alpha = vecdot(r[k - 1], r[k - 1]) / vecdot(s[k], tmp1);

        // x[k] = x[k-1] + alpha * s[k]
        vecscalar(tmp2, alpha, s[k]);
        vecadd(x[k], x[k - 1], tmp2);

        // r[k] = r[k-1] - alpha * A * s[k]
        vecscalar(tmp2, -1 * alpha, tmp1);
        vecadd(r[k], r[k - 1], tmp2);
    }
    veccp(X, x[k]);

    return 0;
}
