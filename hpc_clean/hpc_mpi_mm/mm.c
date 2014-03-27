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
    stripSize = N / numnodes;

    Ablock = matrix_creator(N / numnodes, N);
    Bblock = matrix_creator(N / numnodes, N);
    Cblock = matrix_creator(N / numnodes, N);

    if (myrank == 0) {
        /*start timer*/
        StopwatchStart(watch);

        /* 
         * allocate A, B, and C --- note that you want these to be
         * contiguously allocated.  Workers need less memory allocated.
         */
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
                        tmpB[k] = B[c][r];
                        k++;
                    }
                }
            }
        }

        for (i = 0; i < N; i++) {
            Avett[i] = &tmpA[i * N];
            Bvett[i] = &tmpB[i * N];
        }

        offset = 0;
        numElements = stripSize * N;

        /* send each node its piece of A and B */
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
    /*communicators creation in order to split blocks that will be used in multiplication*/
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
        for (j = 0; j < stripSize * N; j++) {
            Cblock[i][j] = 0.0;
        }
    }

    int rsize, csize;
    double *rbuf;
    double *cbuf;

    MPI_Comm_size(MyComm_row, &rsize);
    MPI_Comm_size(MyComm_col, &csize);

    rbuf = (double *) malloc(rsize * N * sizeof (double));
    cbuf = (double*) malloc(csize * N * sizeof (double));

    MPI_Allgather(Ablock[0], stripSize * N, MPI_DOUBLE, rbuf, stripSize * N, MPI_DOUBLE, MyComm_row);
    MPI_Allgather(Bblock[0], stripSize * N, MPI_DOUBLE, cbuf, stripSize * N, MPI_DOUBLE, MyComm_col);

    /* ripristina la versione a matrice */
    //    double *tmpAA, *tmpBB;
    double **AAblock, **BBblock;

    //    tmpAA = (double *) malloc(rsize * N * sizeof (double));
    //    AAblock = (double **) malloc(sizeof (double *) * rsize);
    //    for (i = 0; i < rsize; i++)
    //        AAblock[i] = &tmpAA[i * N];
    AAblock = matrix_creator(rsize, N);

    //    tmpBB = (double *) malloc(csize * N * sizeof (double));
    //    BBblock = (double **) malloc(sizeof (double *) * csize);
    //    for (i = 0; i < csize; i++)
    //        BBblock[i] = &tmpBB[i * N];
    BBblock = matrix_creator(csize, N);

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

    int l, m;
    for (i = 0; i < N / numnodes; i++) {
        m = 0;
        for (l = 0; l < N / (numnodes / 2); l++) {
            for (j = 0; j < N / (numnodes / 2); j++) {
                for (k = 0; k < N; k++) {
                    Cblock[i][m] += AAblock[l][k] * BBblock[j][k];
                }
                m++;
            }
        }
    }

    /* master receives from workers  -- note could be done via MPI_Gather */
    int gsize;
    double *Carray;
    
    if (myrank == 0) {
        MPI_Comm_size(MPI_COMM_WORLD, &gsize);
        Carray = (double *) malloc(gsize * N * sizeof (double));
    }

    MPI_Gather(Cblock[0], stripSize * N, MPI_DOUBLE, Carray, stripSize * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {

        /* trasform Carray into matrix C */
        k = 0;
        for (i = 0; i < N; i = i + N / (numnodes / 2)) {
            for (j = 0; j < N; j = j + N / (numnodes / 2)) {
                for (r = i; r < i + N / (numnodes / 2); r++) {
                    for (c = j; c < j + N / (numnodes / 2); c++) {
                        C[r][c] = Carray[k];
                        k++;
                    }
                }
            }
        }

        /*stopwatch stop*/
        StopwatchStop(watch);
        StopwatchPrintWithComment("Master total time: %f\n", watch);

        /* print out matrix here, if I'm the master */
        printmatrix(N, N, C);

        freematrix(N, A);
        freematrix(N, B);
        freematrix(N, C);
        free(tmpA);
        free(tmpB);
        free(Avett);
        free(Bvett);
        free(Carray);
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