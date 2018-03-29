#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include "utils/vecmat.h"

#define N 10
#define MAX 100
#define EPS  (1.0e-8)

bool vecIsZero(double vec[])
{
    int i;
    double sum = 0.0;
    for (i = 0; i < N; i++)
        sum += vec[i] * vec[i];
    return (sum < EPS)? true : false;
}

// Preconditioning
void point_jakobi(double A[N][N], double b[])
{
    int i, j;

    double M[N][N];
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            M[i][j] = (i == j)? (double) (1 / A[i][j]) : 0;
        }
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] *= M[i][i];
        }
    }
    veccp(b, matvec(M, b));
}

void cg(double A[N][N], double b[], double** x)
{
    int k = 0; // index
    double *r = malloc(sizeof(double) * N); // 残差ベクトル
    double *s = malloc(sizeof(double) * N); // 最も新しく加えられた基底
    double rr0, rr1; // norm of r (rr0: older, rr1:newer)
    double alpha, beta;

    // r[0] = b - A * x[0] = b
    veccp(r, b);
    while (!vecIsZero(r)) {
        k++;

        if (k == 1) {
            // s[1] = r[0]
            veccp(s, r);
            rr1 = vecdot(r, r);
        } else {
            rr1 = vecdot(r, r);
            beta = rr1 / rr0;
            // sを更新
            // s[k] = r[k-1] + beta * s[k-1]
            s = vecadd(r, vecscalar(beta, s));
        }

        // alpha = (r[k-1] * r[k-1]) / (s[k] * A * s[k])
        double *As = matvec(A, s);
        // printf("rr1 = %f, dot(r, r) = %f\n", rr1, vecdot(r, r));
        // note: 分子はvecdot(r, r)よりrr1の方が収束が早い(桁落ちのため?)
        alpha = (double) (rr1 / vecdot(s, As));

        // xを更新
        // x[k] = x[k-1] + alpha * s[k]
        *x = vecadd(*x, vecscalar(alpha, s));

        // rを更新
        // r[k] = r[k-1] - alpha * A * s[k]
        r = vecsub(r, vecscalar(alpha, As));
        //printVec(*x, k);
        rr0 = rr1;
    }
    free(r); free(s);

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
    double *x = malloc(sizeof(double) * N); // 解
    vecfill(x, 0);

    point_jakobi(A, b);
    cg(A, b, &x);

    int i;
    for (i = 0; i < N; i++) {
        printf("x[%d] = %2g\n", i, x[i]);
    }
    //printVec(x, 0);

    return 0;
}
