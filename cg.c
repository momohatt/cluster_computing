#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <math.h>
#include "utils/vecmat.h"
#include "utils/mat_decomp.h"

#define N 10
#define MAX 100
#define EPS  (1.0e-8)

bool vecIsZero(double vec[])
{
    return (vecnorm(vec) < EPS)? true : false;
}

// Preconditioning
void point_jacobi(double A[N][N], double b[])
{
    int i, j;

    double M[N][N];
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            M[i][j] = (i == j)? (double) (1 / A[i][j]) : 0;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            A[i][j] *= M[i][i];
    matvec(b, M, b);
}

void ssor(double A[N][N], double b[], double omega)
{
    if (!(omega > 0 && omega < 2)) {
        printf("omega should be within range of (0, 2)\n");
        exit(1);
    }
    int i, j, k;
    double d[N], L[N][N]; // Mはいらない

    for (i = 0; i < N; i++) {
        d[i] = A[i][i] - 2.0;
        for (j = 0; j < N; j++) {
            if (i == j)
                L[i][j] = 1.0;
            else
                L[i][j] = (i > j)? A[i][j] : 0;
        }
    }

    vecscalar(d, 1.0 / omega, d); // d = (1 / omega) d
    double DL[N][N];
    double DLt[N][N]; // DL^T
    mat_diag_add(DL, L, d); // DL = (1 / omega) d + L
    double M_rev_A[N][N]; // M^{-1} * A
    double M_rev_b[N]; // M^{-1} * b

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            // d[i]は使えない(1 / omega倍した)のでA[i][i] - 2.0で代用
            DLt[i][j] = DL[j][i] / (A[i][i] - 2.0);
        }
    }

    //printf("D + L = \n");
    //printMat(DL);

    // A, bにM^{-1}を適用していく
    // forward_substitution(DL, [Aのk列目], [Uxのk列目])
    double Ux[N][N]; // L * Ux = A
    for (k = 0; k < N; k++) {
        for (i = 0; i < N; i++) {
            double tmp = A[i][k];
            for (j = 0; j < i; j++) {
                tmp -= DL[i][j] * Ux[j][k];
            }
            Ux[i][k] = tmp / DL[i][i];
        }
    }
    // backward_substitution(DLt, [Uxのk列目], [M_rev_Aのk列目])
    for (k = 0; k < N; k++) {
        for (i = N - 1; i >= 0; i--) {
            double tmp = 0.0;
            for (j = i + 1; j < N; j++) {
                tmp += DLt[i][j] * M_rev_A[j][k];
            }
            M_rev_A[i][k] = (Ux[i][k] - tmp) / DLt[i][i];
        }
    }

    lu_substitution(DL, DLt, b, M_rev_b);

    matcp(A, M_rev_A);
    veccp(b, M_rev_b);
    vecscalar(b, (2.0 - omega) / omega, b);
    matscalar(A, (2.0 - omega) / omega, A);
    printMat(A);
    printVec(b, 0, "b");
    exit(1);
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

void read_input(double A[N][N], double b[])
{
    int i, j;
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            scanf("%lf", &A[i][j]);
    for (i = 0; i < N; i++)
        scanf("%lf", &b[i]);
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

    if (strcmp(argv[1], "jacobi") == 0) {
        printf("CG method with point Jacobi preconditioning\n");
        read_input(A, b);
        point_jacobi(A, b);
    } else if (strcmp(argv[1], "ssor") == 0) {
        if (argc < 3) {
            printf("usage: ./[bin] ssor [omega]\n");
            return 1;
        }
        printf("CG method with SSOR preconditioning\n");
        read_input(A, b);
        ssor(A, b, atof(argv[2]));
    } else if (strcmp(argv[1], "none") == 0) {
        printf("CG method with no preconditioning\n");
        read_input(A, b);
    } else {
        printf("invalid method. usage: ./[bin] [ jacobi | ssor | none ] ([omega])\n");
        return 0;
    }
    cg(A, b, x);

    for (i = 0; i < N; i++) {
        printf("x[%d] = %f\n", i, x[i]);
    }

    return 0;
}
