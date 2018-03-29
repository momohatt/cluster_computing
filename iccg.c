#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "utils/vecmat.h"
#include "utils/mat_decomp.h"

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
void icres(double L[N][N], double d[], double r[], double u[])
{
    int i, j, k;
    double ldlt[N][N]; // L * D * L^T
    double A[N][N + 1]; // ldlt ++ r
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            ldlt[i][j] = L[i][j] * d[i];
        }
    }

    printMat(L);
    transpose(L);
    matmat(ldlt, ldlt, L);
    printMat(ldlt);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] = ldlt[i][j];
        }
        A[i][N] = r[i];
    }

    gaussian_elimination(A, u);

    return;
}

void iccg(double A[N][N], double b[], double *x[])
{
    int i, j, k;
    double *r = malloc(sizeof(double) * N); // 残差ベクトル
    double *s = malloc(sizeof(double) * N); // 方向ベクトル
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

//    double e = 0.0; // error value
    while(!vecIsZero(r)) {
        k++;

        if (k == 1) {
            // s[0] = (LDL^t)^{-1} r[0]
            icres(L, d, r, s);
            rr0 = vecdot(r, s);
        } else {
            icres(L, d, r, r2);
            rr1 = vecdot(r, r2);
            beta = rr1 / rr0;
            s = vecadd(r2, vecscalar(beta, s));
        }
        double *As = matvec(A, s);

        // alpha = r * (LDL^t)^{-1} * r / (s * A * s)
        alpha = rr1 / vecdot(s, As);

        // x, rの更新
        *x = vecadd(*x, vecscalar(alpha, s));
        r = vecsub(r, vecscalar(alpha, As));

//        e = sqrt(rr1);
//        if (e < EPS) {
//            k++;
//            break;
//        }

        // printVec(r, k);
        // rr0, rr1の更新
        rr0 = rr1;
    }
    free(r); free(s); free(r2);
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

//    // Test ICRes
//    double L[N][N];
//    double d[N];
//    ic(A, L, d);
//
//    // L d Lt = A? : OK
//    double r[N] = {1.0, -1.0, 2.0, -2.0, 3.0, -3.0, 4.0, -4.0, 5.0, -5.0};
//    icres(L, d, r, x);

    return 0;
}
