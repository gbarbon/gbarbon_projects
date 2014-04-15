/* 
 * File:   cannon.c
 * Author: asus
 *
 * Created on 14 aprile 2014, 10.56
 */

#include "header.h"
#define TAG 13

int main(int argc, char** argv) {

    double **A, **B, **C, *tmpA, *tmpB, **Ablock, **Bblock;
    double startTime, endTime;
    int nblock, recv, indexR, index, myrank, numnodes, N, i, j, k, r, c;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    N = atoi(argv[1]);
    nblock = numnodes - 1;

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == 0) {
        printf("Printf atoi N: %d\n", N);
        printf("Printf numnodes: %d\n", numnodes);

        A = matrix_creator(N, N);

        B = matrix_creator(N, N);

        C = matrix_creator(N, N);

        printf("Myrank is %d. A,B,C allocated\n", myrank);
    }

    //debug

    if (myrank != 0) {
        A = matrix_creator(N / nblock, N);

        B = matrix_creator(N / nblock, N);

        C = matrix_creator(N / nblock, N);
    }

    if (myrank == 0) {
        // initialize A and B
        simple_matrix_init(A, N);
        simple_matrix_init(B, N);
        
        /* suddivisione in blocchi della matrice */
        tmpA = (double *) malloc(sizeof (double) * N * N);
        tmpB = (double *) malloc(sizeof (double) * N * N);
        Ablock = (double **) malloc(sizeof (double *) * N);
        Bblock = (double **) malloc(sizeof (double *) * N);

        int k = 0;
        for (i = 0; i < N; i = i + N / (numnodes / 2)) {
            for (j = 0; j < N; j = j + N / (numnodes / 2)) {
                for (r = i; r < i + N / (numnodes / 2); r++) {
                    for (c = j; c < j + N / (numnodes / 2); c++) {
                        tmpA[k] = A[r][c];
                        tmpB[k] = B[c][r];
                        k++;
                    }
                }
            }
        }

        for (i = 0; i < N; i++) {
            Ablock[i] = &tmpA[i * N];
            Bblock[i] = &tmpB[i * N];
        }
        
        printmatrix(N, N, Ablock);
        
    }
    else {

    }
}