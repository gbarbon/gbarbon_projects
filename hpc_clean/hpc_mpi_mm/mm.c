/*
 * File:   mm.c
 * Author: 
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
    double **A, **B, **C, **Ablock, **Bblock, **Cblock;
    int numElements, offset, stripSize, myrank, numnodes, N, i, j, k, r, c;

    /*stopwatch*/
    Stopwatch watch = StopwatchCreate();

    /*MPI initialization*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    MPI_Comm MyComm_row;
    MPI_Comm MyComm_col;
    int* coo;

    N = atoi(argv[1]);
    stripSize = N / numnodes;
//    numElements = stripSize * N;
    numElements = (N*N)/numnodes;

    if (myrank == 0) {
        /*start timer*/
        StopwatchStart(watch);

        /*creates matrices*/
        A = matrix_creator(N, N);
        B = matrix_creator(N, N);
        C = matrix_creator(N, N);

        /*initializes A and B randomly*/
        simple_matrix_init(A, N);
        simple_matrix_init(B, N);
        /*printmatrix(N,N,A);*/

        /* transpose B */
        matrix_transposer(N, B);
        /*printmatrix(N,N,B);*/

        /* transform A & B into vectors */
        double * Aarray = matrix_vectorizer(N,N,A);
        double * Barray = matrix_vectorizer(N,N,B);
        
        /* MISSING: split A & B in the correct way */

        /* send each node its piece of A and B */
        offset = 0;
        for (i = 0; i < numnodes; i++) {
            MPI_Send(&Aarray[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(&Barray[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            offset += numElements;
        }
        
        free(Aarray);
        free(Barray);
    }

    Ablock = matrix_creator(stripSize, N);
    Bblock = matrix_creator(stripSize, N);
    Cblock = matrix_creator(stripSize, N);

    /*receive my part of A and B*/
    MPI_Recv(Ablock[0], numElements, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(Bblock[0], numElements, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    /* coords computation */
    coo = coordinate(myrank, numnodes);

    /*communicators creation in order to split blocks that will be used in multiplication*/
    MPI_Comm_split(MPI_COMM_WORLD, coo[0], myrank, &MyComm_row);
    MPI_Comm_split(MPI_COMM_WORLD, coo[1], myrank, &MyComm_col);

    int rsize, csize;
    MPI_Comm_size(MyComm_row, &rsize);
    MPI_Comm_size(MyComm_col, &csize);

    double *rbuf = (double *) malloc(rsize * N * sizeof (double));
    double *cbuf = (double*) malloc(csize * N * sizeof (double));

    MPI_Allgather(Ablock[0], numElements, MPI_DOUBLE, rbuf, numElements, MPI_DOUBLE, MyComm_row);
    MPI_Allgather(Bblock[0], numElements, MPI_DOUBLE, cbuf, numElements, MPI_DOUBLE, MyComm_col);

    /* restore matrix version */
    /*MAY NOT BE CORRECT*/
    double **AAblock, **BBblock;
    AAblock = devectorizer(rsize, N, rbuf);
    BBblock = devectorizer(csize, N, cbuf);
    
    int l, m;
    for (i = 0; i < N / numnodes; i++) {
        m = 0;
        for (l = 0; l < N / (numnodes / 2); l++) {
            for (j = 0; j < N / (numnodes / 2); j++) {
                Cblock[i][m] = 0.0;
                for (k = 0; k < N; k++) {
                    Cblock[i][m] += AAblock[l][k] * BBblock[j][k];
                }
                m++;
            }
        }
    }

    /* master receives from workers  -- note could be done via MPI_Gather */
    double *Carray;

    if (myrank == 0) {
        int gsize;
        MPI_Comm_size(MPI_COMM_WORLD, &gsize);
        Carray = (double *) malloc(gsize * N * sizeof (double));
    }

    MPI_Gather(Cblock[0], numElements, MPI_DOUBLE, Carray, numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        /* trasform Carray into matrix C */
        /*MAY NOT BE CORRECT*/
        C = devectorizer(N, N, Carray);

        /*stopwatch stop*/
        StopwatchStop(watch);
        StopwatchPrintWithComment("Master total time: %f\n\n", watch);

        /* print out matrix here, if I'm the master */
        printmatrix(N, N, C);

        freematrix(N, A);
        freematrix(N, B);
        freematrix(N, C);
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