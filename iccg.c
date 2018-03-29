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

bool vecIsZero(double vec[])
{
    int i;
    double sum = 0.0;
    for (i = 0; i < N; i++)
        sum += vec[i] * vec[i];
    return (sum < EPS)? true : false;
}

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

// u = (LDL^T)^{-1} * r を計算
// (LDL^T) u = r をガウス消去法で解く
void icres(double L[N][N], double d[N], double r[N], double u[N])
{
    int i, j, k;
    double ldlt[N][N]; // A = L * D * L^T
    double A[N][N + 1]; // ldlt ++ r
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            ldlt[i][j] = L[i][j] * d[i];
        }
    }
    transpose(L);
    matmat(ldlt, ldlt, L);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] = ldlt[i][j];
        }
        A[i][N] = r[i];
    }

    // Gauss elimination (forward elimination)
    for (k = 0; k < N - 1; k++) {
        double akk = A[k][k];
        for (i = k + 1; i < N; i++) {
            double aik = A[i][k];
            for (j = k; j <= N; j++) {
                A[i][j] -= aik * (A[k][j] / akk);
            }
        }
    }

    // Gauss elimination (backward substitution)
    u[N - 1] = A[N - 1][N] / A[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {
        double au = 0.0;
        for (int j = i + 1; j < N; j++) {
            au += A[i][j] * u[j];
        }
        u[i] = (A[i][N] - au) / A[i][i];
    }
    return;
}

void iccg(double A[N][N], double b[], double *x[])
{
    int i, j, k;
    double *r = malloc(sizeof(double) * N); // 残差ベクトル
    double *s = malloc(sizeof(double) * N); // 方向ベクトル
    double *y = malloc(sizeof(double) * N); // y = Ap
    double *r2 = malloc(sizeof(double) * N); // r2 = (LDL^T)^{-1} * r

    double rr0, rr1; // rr0 = norm(r[k-1]), rr1 = norm(r[k])
    double alpha, beta;

    double d[N]; // 対角行列Dの対角成分
    double L[N][N];

    matfill(L, 0);
    vecfill(*x, 0);
    ic(A, L, d);

    // r[0] = b - A * x[0] = b[0]
    veccp(r, b);

    // s[0] = (LDL^t)^{-1} r[0]
    icres(L, d, r, s);

    rr0 = vecdot(r, s);

    double e = 0.0; // error value
    for (k = 0; k < MAXITER; k++) {
        // y = A * s
        y = matvec(A, s);

        // alpha = r * (LDL^t)^{-1} * r / (s * A * s)
        alpha = rr0 / vecdot(s, y);

        // x, rの更新
        double tmp[N];
        *x = vecadd(*x, vecscalar(alpha, s));
        r = vecsub(r, vecscalar(alpha, y));
        icres(L, d, r, r2);
        rr1 = vecdot(r, r2);
        e = sqrt(rr1);
        if (e < EPS) {
            k++;
            break;
        }

        beta = rr1 / rr0;
        s = vecadd(r2, vecscalar(beta, s));

        rr0 = rr1;
        printVec(r, k);
    }
    free(r); free(s); free(y); free(r2);
    return;
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
    double *x = malloc(sizeof(double) * N);

    int i, j;

    iccg(A, b, &x);

    for (i = 0; i < N; i++) {
        printf("x[%d] = %2g\n", i, x[i]);
    }

    return 0;
}
