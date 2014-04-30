/*
 * File:   mm.c
 * Author:
 *
 * Created on 12 marzo 2014, 15.20
 */

#include "header.h"
#include "MpiStopwatch.h"
#define TAG 13

/**
 * Calculates the coordinates of each process, in order to let it return
 * the correct rows and cols from other processes
 *
 * @param procNum
 * @param totalProc
 * @return
 */
int* coordinate(int procNum, int totalProc) {
    int* coord = (int*) calloc(2, sizeof (int));
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    return coord;
}

/**
 * Send rows and cols of A and B to the interested worker
 *
 * @param A     Pointer to the first nxn matrix
 * @param B     Pointer to the second nxn matrix
 * @param offset        Dimension of the "local" matrix of the worker
 * @param n     Dimension of the A & B matrices
 */
void master_sender(double** A, double** B, int offset, int n) {
    int i, j, x, y, worker = 0;
    double * tempA, * tempB, ** Ablock, ** Bblock;

    /*stopwatch start*/
    Stopwatch watch = StopwatchCreate();
    StopwatchStart(watch);

    /**/
    Ablock = matrix_creator(offset, offset);
    Bblock = matrix_creator(offset, offset);

    for (i = 0; i < n; i += offset)
        for (j = 0; j < n; j += offset) {

            /**/
            for (x = i; x < offset + i; x++)
                for (y = j; y < offset + j; y++) {
                    Ablock[x - i][y - j] = A[x][y];
                    Bblock[x - i][y - j] = B[x][y];
                }
            /*printmatrix(offset,offset,Ablock);*/

            /*vectorize the two pieces of matrices in order to send them*/
            tempA = matrix_vectorizer(offset, offset, Ablock);
            tempB = matrix_vectorizer(offset, offset, Bblock);
            //printvector(offset*offset,tempA);

            printf("Worker:%d\n", worker);

            /*MPI send of rows and cols. Notice tag 0 for rows, tag 1 for cols*/
            MPI_Send(tempA, offset * offset, MPI_DOUBLE, worker, 0, MPI_COMM_WORLD);
            MPI_Send(tempB, offset * offset, MPI_DOUBLE, worker, 1, MPI_COMM_WORLD);

            printf("Worker:%d\n", worker);

            free(tempA);
            free(tempB);
            worker++; /*increment worker number*/
        }

    /*stopwatch stop*/
    StopwatchStop(watch);
    StopwatchPrintWithComment("Time in master_sender function is: %f\n", watch);
    free(watch);
    free(Ablock);
    free(Bblock);
}

int main(int argc, char *argv[]) {
    int numElements, stripSize, myrank, numnodes, N, i, j, k, x;
    double **A, **B, **C, *Ablock, *Bblock, **Cblock;

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

    //there we may also calculate "lateral" dimension of block:
    int dim = N / sqrt(numnodes);

    /* number of elements for each process*/
    //numElements = (N * N) / numnodes;
    numElements = dim*dim;

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
        matrix_transposer(N, B); /* transpose B */
        printf("\n\n");

        master_sender(A, B, dim, N);
    }

    Ablock = (double *) malloc(numElements * sizeof (double));
    Bblock = (double *) malloc(numElements * sizeof (double));
    Cblock = matrix_creator(dim, dim);

    /*receive my part of A and B*/
    MPI_Recv(Ablock, numElements, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(Bblock, numElements, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    printf("Recv fatte\n");

    /* coords computation */
    coo = coordinate(myrank, numnodes);

    printf("coo[0]: %d coo[1]: %d\n", coo[0], coo[1]);

    /*communicators creation in order to split blocks that will be used in multiplication*/
    MPI_Comm_split(MPI_COMM_WORLD, coo[0], myrank, &MyComm_row);
    //MPI_Comm_split(MPI_COMM_WORLD, myrank, coo[0], &MyComm_row);
    MPI_Comm_split(MPI_COMM_WORLD, coo[1], myrank, &MyComm_col);
    //MPI_Comm_split(MPI_COMM_WORLD, myrank, coo[1], &MyComm_col);

    printf("Split fatto\n");

    /*rsize is the number of square blocks in a row block*/
    /*csize is the number of square blocks in a column block*/
    int rsize, csize;
    MPI_Comm_size(MyComm_row, &rsize);
    MPI_Comm_size(MyComm_col, &csize);

    printf("Myrank %d, Comm_size fatta. rsize= %d\n", myrank, rsize);
    printf("Myrank %d, Comm_size fatta. csize= %d\n", myrank, csize);

    /*rsize times numElements is the number of elements of a row block*/
    double *rbuf = (double *) malloc(rsize * numElements * sizeof (double));
    double *cbuf = (double *) malloc(csize * numElements * sizeof (double));

    printf("rubuf= %d\n", rsize * numElements);

    MPI_Allgather(Ablock, numElements, MPI_DOUBLE, rbuf, numElements, MPI_DOUBLE, MyComm_row);
    MPI_Allgather(Bblock, numElements, MPI_DOUBLE, cbuf, numElements, MPI_DOUBLE, MyComm_col);

    printf("AllGather fatta\n");

    /* restore matrix version */
    double ** BBtest = matrix_creator(N, dim);
    x = 0;
    for (i = 0; i < N; i++)
        for (j = 0; j < dim; j++) {
            BBtest[i][j] = cbuf[x];
            x++;
        }

    double ** AAtest = matrix_creator(dim, N);
    x = 0;
    for (k = 0; k < N; k += dim)
        for (i = 0; i < dim; i++)
            for (j = k; j < k + dim; j++) {
                AAtest[i][j] = rbuf[x];
                x++;
            }

    int l, m;
    for (l = 0; l < dim; l++)
        for (j = 0; j < dim; j++) {
            Cblock[l][j] = 0.0;
            for (k = 0; k < N; k++) {
                Cblock[l][j] += AAtest[l][k] * BBtest[k][j];
            }
        }

    double * cvector = matrix_vectorizer(dim, dim, Cblock);

    /* master receives from workers  -- note could be done via MPI_Gather */
    double *Carray;

    if (myrank == 0) {
        int gsize;
        MPI_Comm_size(MPI_COMM_WORLD, &gsize);
        Carray = (double *) malloc(N * N * sizeof (double));
    }

    MPI_Gather(cvector, numElements, MPI_DOUBLE, Carray, numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("Fino a Gather OK\n");

    if (myrank == 0) {
        /* trasform Carray into matrix C */

        double ** C = matrix_creator(N, N);
        int i, j, x, y, el = 0;
        /*base point (in final matrix) scrolling*/
        /*i.e. 0,0 - 0,2 - 2,0 - 2,2 with N=4 and nproc=4 */
        for (i = 0; i < N; i += dim)
            for (j = 0; j < N; j += dim)
                /*block scorlling, dim: block dimension*/
                for (x = i; x < dim + i; x++)
                    for (y = j; y < dim + j; y++) {
                        C[x][y] = Carray[el];
                        el++;
                    }

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

    free(cvector);
    free(Ablock);
    free(Bblock);
    freematrix(N / numnodes, Cblock);
    free(coo);
    freematrix(dim, AAtest);
    freematrix(N, BBtest);
    free(rbuf);
    free(cbuf);
    free(watch);

    MPI_Finalize();
    return 0;
}