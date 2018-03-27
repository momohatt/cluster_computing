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

int main (int argc, char *argv[])
{
    int k = 0; // index
    double x[N]; // 解(暫定)
    double r[N]; // 残差ベクトル
    double s[N]; // 最も新しく加えられた基底
    double r_norm_cache[MAX];

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
    vecfill(x, 0);

    int i, j;

    // preconditioning
    double M[N][N];
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            M[i][j] = (i == j)? (double) (1 / A[i][j]) : 0;
        }
    }
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] = M[i][i] * A[i][j];
        }
    }
    matvec(b, M, b);


    veccp(r, b);
    while (!vecIsZero(r)) {
        k++;
        double tmp1[N], tmp2[N];

        if (k == 1) {
            // r[k] = b - A * x[k]
//            matvec(tmp1, A, x);
//            vecscalar(tmp1, -1, tmp1);
//            vecadd(r, b, tmp1);
//
            // s[1] = r[0]
            veccp(s, r);
            r_norm_cache[0] = vecdot(r, r);
        } else {
            r_norm_cache[k - 1] = vecdot(r, r);
            double beta = (double) r_norm_cache[k - 1] / r_norm_cache[k - 2];
            // sを更新
            // s[k] = r[k-1] + beta * s[k-1]
            vecscalar(tmp1, beta, s);
            vecadd(s, r, tmp1);
        }

        // alpha = (r[k-1] * r[k-1]) / (s[k] * A * s[k])
        matvec(tmp1, A, s);
        double alpha = (double) (vecdot(r, r) / vecdot(s, tmp1));

        // xを更新
        // x[k] = x[k-1] + alpha * s[k]
        vecscalar(tmp2, alpha, s);
        vecadd(x, x, tmp2);

        // rを更新
        // r[k] = r[k-1] - alpha * A * s[k]
        vecscalar(tmp2, -1 * alpha, tmp1);
        vecadd(r, r, tmp2);
        printVec(x, k);
    }
    for (i = 0; i < N; i++) {
        printf("x[%d] = %2g\n", i, x[i]);
    }
    //printVec(x, 0);

    return 0;
}
