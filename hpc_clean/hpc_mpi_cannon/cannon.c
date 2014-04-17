/* 
 * File:   cannon.c
 * Author: asus
 *
 * Created on 14 aprile 2014, 10.56
 */

#include "header.h"
#define TAG 13

int** coordinate(int totalProc) {
    int i, var;

    int **coord = (int **) malloc(totalProc * sizeof (int*));
    for (i = 0; i < totalProc; i++)
        coord[i] = (int*) calloc(2, sizeof (int));

    var = 3;

    for (i = 0; i < totalProc; i++) {
        coord[i][0] = i / var;
        coord[i][1] = i % var;
    }

    return coord;
}

void skewing_row(double ** M, int n) {
    int i, j, k, index;

    double **r_swap = (double **) malloc(sizeof (double*) * n);
    for (i = 0; i < n; i++) {
        r_swap[i] = M[i];
    }

    k = 1;
    for (i = sqrt(n); i < n; i = i + sqrt(n)) {
        for (j = 0; j < sqrt(n); j++) {
            index = (j + k) % (int) sqrt(n);
            M[i + j] = r_swap[i + index];
        }
        k++;
    }

    //freematrix(n, r_swap);
}

void skewing_column(double ** M, int n) {
    int i, j, k, index;

    double **c_swap = (double **) malloc(sizeof (double*) * n);
    for (i = 0; i < n; i++) {
        c_swap[i] = M[i];
    }

    k = 1;
    for (i = 1; i < n; i = i + sqrt(n)) {
        for (j = 0; j < sqrt(n); j++) {
            index = (j + k) % (int) sqrt(n);
            M[i + (int) sqrt(n) * j] = c_swap[i + (int) sqrt(n) * index];
        }
        k++;
    }

    //freematrix(n, c_swap);
}

double ** matrix_block(double ** matrix, int block, int n) {
    double *tmpM, **Mblock;
    int i, j, r, c, k;

    tmpM = (double *) malloc(sizeof (double) * n * n);
    Mblock = (double **) malloc(sizeof (double *) * n);

    //Mblock = matrix_creator(n, n);

    k = 0;
    for (i = 0; i < n; i = i + n / (block / 2)) {
        for (j = 0; j < n; j = j + n / (block / 2)) {
            for (r = i; r < i + n / (block / 2); r++) {
                for (c = j; c < j + n / (block / 2); c++) {
                    tmpM[k] = matrix[r][c];
                    k++;
                }
            }
        }
    }

    for (i = 0; i < n; i++) {
        Mblock[i] = &tmpM[i * n];
    }

    return Mblock;
}

int main(int argc, char** argv) {

    double **A, **B, **C, *tmpA, *tmpB, **Ablock, **Bblock;
    double startTime, endTime;
    int **coo, nblock, recv, indexR, index, myrank, numnodes, N, i, j, k, r, c;

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

        // suddivisione in blocchi della matrice
        Ablock = (double **) malloc(sizeof (double *) * N);
        Bblock = (double **) malloc(sizeof (double *) * N);

        Ablock = matrix_block(A, nblock, N);
        Bblock = matrix_block(B, nblock, N);

        printmatrix(N, N, Ablock);
        printmatrix(N, N, Bblock);

        // skewing
        skewing_row(Ablock, N);
        skewing_column(Bblock, N);

        printmatrix(N, N, Ablock);
        printmatrix(N, N, Bblock);

    }



    MPI_Finalize();

    return 0;
}