#include "mat_decomp.h"

void forward_substitution(const double L[N][N], const double b[], double x[])
{
    for (int i = 0; i < N; i++) {
        double tmp = b[i];
        for (int j = 0; j < i; j++) {
            tmp -= L[i][j] * x[j];
        }
        x[i] = tmp / L[i][i];
    }
}

void backward_substitution(const double U[N][N], const double b[], double x[])
{
    for (int i = N - 1; i >= 0; i--) {
        double tmp = 0.0;
        for (int j = i + 1; j < N; j++) {
            tmp += U[i][j] * x[j];
        }
        x[i] = (b[i] - tmp) / U[i][i];
    }
}

void lu_substitution(const double L[N][N], const double U[N][N], const double r[N], double u[N])
{
    double y[N];

    // L * y = rとなるyを求める
    forward_substitution(L, r, y);

    // y = U * uなるuを求める
    backward_substitution(U, y, u);
}

void gaussian_elimination(double A[N][N], double b[], double x[])
{
    // Forward elimination
    // - 対角要素をのぞいた左下要素を全て0にする(上三角行列にする)
    for (int k = 0; k < N - 1; k++) {
        double akk = A[k][k];
        for (int i = k + 1; i < N; i++) {
            double aik = A[i][k];
            for (int j = k; j < N; j++) {
                A[i][j] -= aik * (A[k][j] / akk);
            }
            b[i] -= aik * (b[k] / akk);
        }
    }

    // Back substitution
    x[N - 1] = b[N - 1] / A[N - 1][N - 1];
    for (int i = N - 2; i >= 0; i--) {
        double ax = 0.0;
        for (int j = i + 1; j < N; j++) {
            ax += A[i][j] * x[j];
        }
        x[i] = (b[i] - ax) / A[i][i];
    }
}

void lu(double A[N][N], double L[N][N], double U[N][N])
{
    for (int p = 0; p < N; p++) {
        L[p][p] = A[p][p];
        for (int i = p + 1; i < N; i++) {
            L[i][p] = A[i][p];
            L[p][i] = 0;
        }
        U[p][p] = 1;
        for (int i = p + 1; i < N; i++) {
            U[p][i] = A[p][i] / A[p][p];
            U[i][p] = 0;
        }
        for (int i = p + 1; i < N; i++) {
            for (int j = p + 1; j < N; j++) {
                A[i][j] -= U[p][j] * L[i][p];
            }
        }
    }
}
