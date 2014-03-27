/*
 * File:   mpi_mm.c
 * Author: jian
 *
 * Created on 12 marzo 2014, 15.20
 */

#include "header.h"
#include "MpiStopwatch.h"
#define TAG 13

int* coordinate(int procNum, int totalProc) {
    int* coord = (int*) calloc(2, sizeof (int)); 
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    return coord;
}

int main(int argc, char *argv[]) {
    double **A, **B, **C, *tmpA, *tmpB, **Avett, **Bvett;
    double **Ablock, **Bblock, **Cblock;
    int numElements, offset, stripSize, myrank, numnodes, N, i, j, k, r, c;

    /*stopwatch*/
    Stopwatch watch = StopwatchCreate();    

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    MPI_Comm MyComm_row;
    MPI_Comm MyComm_col;
    int* coo;

    N = atoi(argv[1]);
    numnodes = 4;

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    Ablock = matrix_creator(N / numnodes, N);
    Bblock = matrix_creator(N / numnodes, N);
    Cblock = matrix_creator(N / numnodes, N);

    if (myrank == 0) {

        A = matrix_creator(N, N);
        B = matrix_creator(N, N);
        C = matrix_creator(N, N);
        
        matrix_init(A, N);
        matrix_init(B, N);

        /* suddivisione in blocchi della matrice */
        tmpA = (double *) malloc(sizeof (double) * N * N);
        tmpB = (double *) malloc(sizeof (double) * N * N);
        Avett = (double **) malloc(sizeof (double *) * N);
        Bvett = (double **) malloc(sizeof (double *) * N);

        int k = 0;
        for (i = 0; i < N; i = i + N / (numnodes / 2)) {
            for (j = 0; j < N; j = j + N / (numnodes / 2)) {
                for (r = i; r < i + N / (numnodes / 2); r++) {
                    for (c = j; c < j + N / (numnodes / 2); c++) {
                        tmpA[k] = A[r][c];
                        tmpB[k] = B[r][c];
                        k++;
                    }
                }
            }
        }

        for (i = 0; i < N; i++) {
            Avett[i] = &tmpA[i * N];
            Bvett[i] = &tmpB[i * N];
        }
    }

    stripSize = N / numnodes;

    /* send each node its piece of A and B */
    if (myrank == 0) {
        /*start timer*/
        StopwatchStart(watch);
        
        offset = 0;
        numElements = stripSize * N;

        for (i = 0; i < numnodes; i++) {
            MPI_Send(Avett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(Bvett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            offset += stripSize;
        }
    }
    /*receive my part of A and B*/
    MPI_Recv(Ablock[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(Bblock[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* coords computation */
    coo = coordinate(myrank, numnodes);

    /*MPI_Barrier(MPI_COMM_WORLD);*/
    /*creazione communicatori per la condivisione dei blocchi necessari alla moltiplicazione*/
    MPI_Comm_split(MPI_COMM_WORLD, coo[0], myrank, &MyComm_row);
    MPI_Comm_split(MPI_COMM_WORLD, coo[1], myrank, &MyComm_col);

    /* Let each process initialize C to zero */
    if (myrank == 0) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                C[i][j] = 0.0;
            }
        }
    }
    for (i = 0; i < N / numnodes; i++) {
        for (j = 0; j < N / numnodes; j++) {
            Cblock[i][j] = 0.0;
        }
    }

    int rsize, csize;
    double *rbuf;
    double *cbuf;

    MPI_Comm_size(MyComm_row, &rsize);
    MPI_Comm_size(MyComm_col, &csize);

    rbuf = (double *) malloc(rsize * 4 * sizeof (double));
    cbuf = (double*) malloc(csize * 4 * sizeof (double));

    MPI_Allgather(Ablock[0], stripSize * N, MPI_DOUBLE, rbuf, stripSize * N, MPI_DOUBLE, MyComm_row);
    MPI_Allgather(Bblock[0], stripSize * N, MPI_DOUBLE, cbuf, stripSize * N, MPI_DOUBLE, MyComm_col);

    /* ripristina la versione a matrice */
//    double *tmpAA, *tmpBB;
    double **AAblock, **BBblock;

    //    tmpAA = (double *) malloc(rsize * 4 * sizeof (double));
    //    AAblock = (double **) malloc(sizeof (double *) * rsize);
    //    for (i = 0; i < rsize; i++)
    //        AAblock[i] = &tmpAA[i * N];
    AAblock = matrix_creator(rsize, 4);

    //    tmpBB = (double *) malloc(csize * 4 * sizeof (double));
    //    BBblock = (double **) malloc(sizeof (double *) * csize);
    //    for (i = 0; i < csize; i++)
    //        BBblock[i] = &tmpBB[i * N];
    BBblock = matrix_creator(csize, 4);

    k = 0;
    for (j = 0; j < N; j = j + N / (numnodes / 2)) {
        for (r = 0; r < N / (numnodes / 2); r++) {
            for (c = j; c < j + N / (numnodes / 2); c++) {
                AAblock[r][c] = rbuf[k];
                BBblock[r][c] = cbuf[k];
                k++;
            }
        }
    }

    if (myrank == 0) {
        freematrix(N, A);
        freematrix(N, B);
        freematrix(N, C);
        free(tmpA);
        free(tmpB);
        free(Avett);
        free(Bvett);

        /*stopwatch stop*/
        StopwatchStop(watch);
        StopwatchPrintWithComment("Master total time: %f\n", watch);
    }

    freematrix(N / numnodes, Ablock);
    freematrix(N / numnodes, Bblock);
    freematrix(N / numnodes, Cblock);
    free(coo);
    freematrix(rsize, AAblock);
    freematrix(csize, BBblock);
    free(rbuf);
    free(cbuf);
    free(watch);

    MPI_Finalize();
    return 0;
}