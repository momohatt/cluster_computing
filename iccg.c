#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "utils/vecmat.h"

#define N 10
#define MAXITER 100
#define EPS  (1.0e-8)

// Incomplete Cholesky Decomposition
void ic(double A[N][N], double L[N][N], double d[N])
{
    L[0][0] = A[0][0];
    d[0] = 1.0 / L[0][0];

    int i, j, k;
    for (i = 1; i < N; i++) {
        for (j = 0; j <= i; j++) {
            if (fabs(A[i][j]) < EPS) continue;

            double lld = A[i][j];
            for (k = 0; k < j; k++) {
                lld -= L[i][k] * L[j][k] * d[k];
            }
            L[i][j] = lld;
        }
        d[i] = 1.0 / L[i][i];
    }
}

void icres(const double L[N][N], const double d[N], const double r[N], double u[N])
{
    int i, j;

    // L * y = rとなるyを求める
    double y[N];
    for (i = 0; i < N; i++) {
        double tmp = r[i];
        for (j = 0; j < i; j++) {
            tmp -= L[i][j] * y[j];
        }
        y[i] = tmp / L[i][i];
    }

    // y = (D Lt) uなるuを求める
    for (i = N - 1; i >= 0; i--) {
        double tmp = 0.0;
        for (j = i + 1; j < N; j++) {
            tmp += L[j][i] * u[j];
        }
        u[i] = y[i] - d[i] * tmp;
    }
}

void iccg(double A[N][N], double b[], double x[])
{
    int i, j, k = 0;
    double r[N]; // 残差ベクトル
    double s[N]; // 方向ベクトル
    double r2[N]; // r2 = (LDL^T)^{-1} * r

    double rr0, rr1; // rr0 = norm(r[k-1]), rr1 = norm(r[k])
    double alpha, beta;

    double d[N]; // 対角行列Dの対角成分
    double L[N][N];

    matfill(L, 0);
    vecfill(x, 0);
    ic(A, L, d);

    // r[0] = b - A * x[0] = b[0]
    veccp(r, b);

    double err = 1.0; // 初期値は適当(EPS以上ならなんでもよい)
    while(err > EPS) {
        double tmp[N];
        printVec(x, k, "x");
        k++;

        if (k == 1) {
            // s[0] = (LDL^t)^{-1} r[0]
            icres(L, d, r, s);
            rr1 = vecdot(r, s);
        } else {
            // r2 = (LDL^t)^{-1} r
            icres(L, d, r, r2);
            printVec(r, k, "r");
            //printVec(r2, k, "r2");
            rr1 = vecdot(r, r2);
            beta = (double) (rr1 / rr0);
            vecscalar(tmp, beta, s);
            vecadd(s, r2, tmp);
        }
        double As[N];
        matvec(As, A, s);

        // alpha = r * (LDL^t)^{-1} * r / (s * A * s)
        alpha = (double) (rr1 / vecdot(s, As));

        //printf("alpha = %f\n", alpha);
        //printf("rr0 = %f\n", rr0);
        //printf("rr1 = %f\n", rr1);

        // x, rの更新
        vecscalar(tmp, alpha, s);
        vecadd(x, x, tmp);
        vecscalar(tmp, alpha, As);
        vecsub(r, r, tmp);

        err = sqrt(rr1);

        // rr0, rr1の更新
        rr0 = rr1;
    }
    return;
}

int main (int argc, char *argv[])
{
    double A[N][N];
    double b[N];
    double x[N]; // 解
    vecfill(x, 0);
    int i, j;

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            scanf("%lf", &A[i][j]);
    for (i = 0; i < N; i++)
        scanf("%lf", &b[i]);

    printf("ICCG method\n");
    iccg(A, b, x);

    for (i = 0; i < N; i++) {
        printf("x[%d] = %f\n", i, x[i]);
    }

    return 0;
}
