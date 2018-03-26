#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define TIMES 1000000000

int main(int argc, char **argv)
{
    //long long int niter = 0; // number of iteration
    long long int count = 0, mycount = 0;
    double x, y;
    double pi;
    int myrank, numprocs, proc;
    int master = 0;
    int tag = 123;
    clock_t start_comp, end_comp;
    clock_t start_comm, end_comm;
    if (myrank == master)
        start_comm = clock();
    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    srand((unsigned)time(NULL) + myrank);

//    if (argc <= 1) {
//        fprintf(stderr, "input number_of_iterations\n");
//        MPI_Finalize();
//        exit(-1);
//    }

    //sscanf(argv[1], "%llu", &niter);
    if (myrank == master)
        start_comp = clock();
    long long int i;
    for (i = 0; i < TIMES / numprocs; i++) {
        //srand((unsigned) (time(NULL) + myrank + i));
        x = (double) rand() / RAND_MAX;
        //srand((unsigned) (time(NULL) + myrank + i));
        y = (double) rand() / RAND_MAX;
        if (x * x + y * y <= 1)
            mycount++;
    }
    if (myrank == master)
        end_comp = clock();

    if (myrank == master) {
        count = mycount;
        for (proc = 1; proc < numprocs; proc++) {
            MPI_Recv(&mycount,  // 受信するデータの先頭アドレス
                     1,         // データの個数
                     MPI_REAL,  // データ型
                     proc,      // 送信元プロセスのランク
                     tag,       // メッセージ識別番号
                     MPI_COMM_WORLD,
                     &status);  // 送信プロセスの情報
            count += mycount;
        }
        pi = (double) count / (TIMES / numprocs * numprocs) * 4;
        printf("number of processors = %d, number of trials = %d, pi = %.10f, ", numprocs, TIMES / numprocs * numprocs, pi);
        printf("time of computation: %f sec, ", (double) (end_comp - start_comp) / CLOCKS_PER_SEC);
    } else {
        //printf("Processor %d sending results = %d to master process\n",
        //       myrank, mycount);
        MPI_Send(&mycount,  // 送信するデータの先頭アドレス
                 1,         // データの個数
                 MPI_REAL,  // データ型
                 master,    // 送信先プロセスのランク
                 tag,       // メッセージ識別番号
                 MPI_COMM_WORLD);
    }
    if (myrank == master) {
        end_comm = clock();
        printf("time of communication: %f sec\n", (double) ((end_comm - start_comm) - (end_comp - start_comp)) / CLOCKS_PER_SEC);
    }

    MPI_Finalize();
}
