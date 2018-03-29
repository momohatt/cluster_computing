#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "utils/vecmat.h"

#define N 10
#define MAX 100
#define EPS  (1.0e-8)

bool vecIsZero(double vec[])
{
    return (vecnorm(vec) < EPS)? true : false;
}

void apply_precond(double A[N][N], double b[], double M[N][N])
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            A[i][j] *= M[i][i];
    matvec(b, M, b);
}

// Preconditioning
void point_jacobi(double A[N][N], double b[])
{
    int i, j;

    double M[N][N];
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            M[i][j] = (i == j)? (double) (1 / A[i][j]) : 0;
    apply_precond(A, b, M);
}

void ssor(double A[N][N], double b[], double omega)
{
    if (!(omega > 0 && omega < 2)) {
        printf("omega should be within range of (0, 2)\n");
        return;
    }
    int i, j;
    double M[N][N], d[N], L[N][N];

    for (i = 0; i < N; i++) {
        d[i] = A[i][i];
        for (j = 0; j < N; j++) {
            L[i][j] = (i > j)? A[i][j] : 0;
        }
    }

    vecscalar(d, 1.0 / omega, d);
    double DL[N][N];
    double DLt[N][N];
    mat_diag_add(DL, L, d);

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            DLt[i][j] = DL[j][i];
            DL[i][j] /= d[j];
        }
    }
//    printf("(D + L) * D^{-1} = \n");
//    printMat(DL);
//    printf("(D + L)^T = \n");
//    printMat(DLt);
    matmat(DL, DL, DLt);
    matscalar(M, 1.0 / (2.0 - omega), DL);

    apply_precond(A, b, M);
    return;
}

void cg(double A[N][N], double b[], double x[])
{
    int k = 0; // index
    double r[N]; // 残差ベクトル
    double s[N]; // 最も新しく加えられた基底
    double rr0, rr1; // norm of r (rr0: older, rr1:newer)
    double alpha, beta;

    // r[0] = b - A * x[0] = b
    veccp(r, b);
    while (!vecIsZero(r)) {
        k++;
        double tmp1[N], tmp2[N];

        if (k == 1) {
            // s[1] = r[0]
            veccp(s, r);
            rr1 = vecdot(r, r);
        } else {
            rr1 = vecdot(r, r);
            beta = rr1 / rr0;
            // sを更新
            // s[k] = r[k-1] + beta * s[k-1]
            vecscalar(tmp1, beta, s);
            vecadd(s, r, tmp1);
        }

        // alpha = (r[k-1] * r[k-1]) / (s[k] * A * s[k])
        matvec(tmp1, A, s);
        // printf("rr1 = %f, dot(r, r) = %f\n", rr1, vecdot(r, r));
        // note: 分子はvecdot(r, r)よりrr1の方が収束が早い(桁落ちのため?)
        alpha = (double) (rr1 / vecdot(s, tmp1));

        // xを更新
        // x[k] = x[k-1] + alpha * s[k]
        vecscalar(tmp2, alpha, s);
        vecadd(x, x, tmp2);

        // rを更新
        // r[k] = r[k-1] - alpha * A * s[k]
        vecscalar(tmp2, alpha, tmp1);
        vecsub(r, r, tmp2);
        printVec(x, k, "x");
        rr0 = rr1;
    }
}

int main (int argc, char *argv[])
{
    double A[N][N];
    double b[N];
    double x[N]; // 解
    vecfill(x, 0);
    int i, j;

    if (argc < 2) {
        printf("usage: ./[bin] [ jacobi | ssor | none ] ([omega])\n");
        return 1;
    }

    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            scanf("%lf", &A[i][j]);
    for (i = 0; i < N; i++)
        scanf("%lf", &b[i]);

    if (strcmp(argv[1], "jacobi") == 0) {
        point_jacobi(A, b);
    } else if (strcmp(argv[1], "ssor") == 0) {
        if (argc < 3) {
            printf("usage: ./[bin] ssor [omega]\n");
            return 1;
        }
        ssor(A, b, atof(argv[2]));
    }
    cg(A, b, x);

    for (i = 0; i < N; i++) {
        printf("x[%d] = %f\n", i, x[i]);
    }

    return 0;
}
