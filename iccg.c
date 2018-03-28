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

    for (k = 0; k < N - 1; k++) {
        double akk = A[k][k];
        for (i = k + 1; i < N; i++) {
            double aik = A[i][k];
            for (j = k; j <= N; j++) {
                A[i][j] -= aik * (A[k][j] / akk);
            }
        }
    }
    u[N - 1] = A[N - 1][N] / A[N - 1][N - 1];
    for (i = N - 2; i >= 0; i--) {
        double ax = 0.0;
        for (int j = i + 1; j < N; j++) {
            ax += A[i][j] * u[j];
        }
        u[i] = (A[i][N] - ax) / A[i][i];
    }
    return;
}

void iccg(double A[N][N], double b[N], double x[N])
{
    int i, j, k;
    double r[N]; // 残差ベクトル
    double s[N]; // 方向ベクトル
    double y[N]; // y = Ap
    double r2[N]; // r2 = (LDL^T)^{-1} * r

    double rr0, rr1; // rr0 = norm(r[k-1]), rr1 = norm(r[k])
    double alpha, beta;

    double d[N]; // 対角行列D
    double L[N][N];

    matfill(L, 0);
    vecfill(x, 0);
    ic(A, L, d);

    // r = b - A * x = b
    veccp(r, b);

    // p_0 = (LDL^t)^{-1} r_0
    icres(L, d, r, s);

    rr0 = vecdot(r, s);

    double e = 0.0; // error?
    for (k = 0; k < MAXITER; k++) {
        // y = A * s
        matvec(y, A, s);

        // alpha = r * (LDL^t)^{-1} * r / (s * A * s)
        alpha = rr0 / vecdot(s, y);

        // x, rの更新
        double tmp[N];
        vecscalar(tmp, alpha, s);
        vecadd(x, x, tmp);
        vecscalar(tmp, alpha, y);
        vecsub(r, r, tmp);
        icres(L, d, r, r2);
        rr1 = vecdot(r, r2);
        e = sqrt(rr1);
        if (e < EPS) {
            k++;
            break;
        }

        beta = rr1 / rr0;
        vecscalar(tmp, beta, s);
        vecadd(s, r2, tmp);

        rr0 = rr1;
    }
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
    double x[N];

    int i, j;

    iccg(A, b, x);

    for (i = 0; i < N; i++) {
        printf("x[%d] = %2g\n", i, x[i]);
    }

    return 0;
}
