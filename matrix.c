#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define SIZE 1000
#define TIMES 10000

int main (int argc, char *argv[])
{
    long long int i, j;
    long long int t;
    int a[SIZE][SIZE], b[SIZE], c[SIZE];
    clock_t start, end;
    for (i = 0; i < SIZE; i++) {
        b[i] = 1;
        for (j = 0; j < SIZE; j++) {
            a[i][j] = (i == j)? 1 : 0;
        }
    }

    printf("hoge\n");
    start = clock();
    for (t = 0; t < TIMES; t++) {
        for (i = 0; i < SIZE; i++) {
            for (j = 0; j < SIZE; j++) {
                c[i] += + a[i][j] * b[j];
            }
        }
    }
    end = clock();
    printf("IJ: %f sec\n", (double)((end - start) / CLOCKS_PER_SEC));

    start = clock();
    for (t = 0; t < TIMES; t++) {
        for (j = 0; j < SIZE; j++) {
            for (i = 0; i < SIZE; i++) {
                c[i] += + a[i][j] * b[j];
            }
        }
    }
    end = clock();
    printf("JI: %f sec\n", (double)((end - start) / CLOCKS_PER_SEC));

    int is;
    start = clock();
    for (t = 0; t < TIMES; t++) {
        for (is = 0; is < SIZE; is += 10) {
            for (j = 0; j < SIZE; j++) {
                for (i = is; i < is + 9; i++) {
                    c[i] += + a[i][j] * b[j];
                }
            }
        }
    }
    end = clock();
    printf("IsJI: %f sec\n", (double)((end - start) / CLOCKS_PER_SEC));

    return 0;
}
