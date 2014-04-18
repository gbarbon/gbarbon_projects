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
        coord[i] = (int*) calloc(3, sizeof (int));

    var = sqrt(totalProc);

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
    for (i = 1; i < n - 1; i = i + sqrt(n)) {
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

// determinazione del rank del processo a cui devo inviare il sottoblocco di A posseduto

int getRankRowDest(int rank, int np) {
    int rank_dest;
    int proc_lato = (int) sqrt(np);
    // se siamo ad inizio riga
    if ((rank % proc_lato) == 0)
        rank_dest = rank + (proc_lato - 1);
    else
        rank_dest = rank - 1;

    return rank_dest;
}

// determinazione del rank del processo da cui devo ricevere il sottoblocco di A

int getRankRowMit(int rank, int np) {
    int rank_mit;
    int proc_lato = (int) sqrt(np);
    // se siamo sull'ultima colonna
    if ((rank + 1) % proc_lato == 0)
        rank_mit = rank - (proc_lato - 1);
    else
        rank_mit = rank + 1;

    return rank_mit;
}

// determinazione del rank del processo a cui devo inviare il sottoblocco di B

int getRankColDest(int rank, int np) {
    int rank_dest;
    int proc_lato = (int) sqrt(np);

    if (rank < proc_lato)
        rank_dest = np - (proc_lato - rank);
    else
        rank_dest = rank - proc_lato;

    return rank_dest;
}

// determinazione del rank del processo da cui ricevere il sottoblocco di B

int getRankColMit(int rank, int np) {
    int rank_mit;
    int proc_lato = (int) sqrt(np);

    if (rank >= (np - proc_lato))
        rank_mit = (rank + proc_lato) % np;
    else
        rank_mit = rank + proc_lato;

    return rank_mit;
}

int main(int argc, char** argv) {

    double **A, **B, **C, *tmpA, *tmpB, **Ablock, **Bblock;
    double startTime, endTime;
    int **coo, nblock, stripSize, numElements, offset, myrank, numnodes, N, i, j, k;

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
    } else {
        A = matrix_creator(N / nblock, N);

        B = matrix_creator(N / nblock, N);

        C = matrix_creator(N / nblock, N);
    }

    stripSize = N / nblock;

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

        // send
        offset = 0;
        numElements = stripSize * N;

        for (i = 1; i <= nblock; i++) {
            MPI_Send(Ablock[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(Bblock[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);

            offset += stripSize;
        }

        printf("Myrank is %d. Must be 0. Pieces of A and B sent.\n", myrank);

    } else {
        // receive my part of A and B
        numElements = stripSize * N;

        MPI_Recv(A[0], numElements, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(B[0], numElements, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Myrank is %d. Must NOT be 0. Pieces of A and B received.\n", myrank);

        printmatrix(N / nblock, N, A);
        printmatrix(N / nblock, N, B);
    }

    MPI_Finalize();

    return 0;
}