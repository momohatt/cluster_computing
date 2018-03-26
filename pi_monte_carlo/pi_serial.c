#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define TIMES 1000000000

int main (int argc, char *argv[])
{
    long long int i, cnt = 0;
    clock_t start, end;
    start = clock();
    for (i = 0; i < TIMES; i++) {
        //if (i % (TIMES / 100) == 0)
        //    printf("%llu: pi = %.10f\n", i, (double) (4 * cnt) / i);
        double px = (double) (rand() + 1) / (RAND_MAX + 2);
        double py = (double) (rand() + 1) / (RAND_MAX + 2);
        if ((px * px + py * py) <= 1)
            cnt++;
    }
    end = clock();
    printf("serial computation\n");
    printf("number of trials = %d, pi = %.10f\n", TIMES, (double) (4 * cnt) / TIMES);
    printf("time: %f sec\n", (double) (end - start) / CLOCKS_PER_SEC);
    return 0;
}
