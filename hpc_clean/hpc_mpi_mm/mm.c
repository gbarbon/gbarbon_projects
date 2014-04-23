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
            /*printvector(offset*offset,tempA);*/

            /*MPI send of rows and cols. Notice tag 0 for rows, tag 1 for cols*/
            MPI_Send(tempA, offset * offset, MPI_DOUBLE, worker, 0, MPI_COMM_WORLD);
            MPI_Send(tempB, offset * offset, MPI_DOUBLE, worker, 1, MPI_COMM_WORLD);

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
    //    double **A, **B, **C, **Ablock, **Bblock, **Cblock;
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
    stripSize = N / numnodes;
    //    numElements = stripSize * N;
    /* number of elements for each process*/
    numElements = (N * N) / numnodes;
    //there we may also calculate "lateral" dimension of block:
    int dim = N / sqrt(numnodes);

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
        printf("Matrix A:\n");
        printmatrix(N, N, A);
        printf("\n\nMatrix B:\n");
        printmatrix(N, N, B);
        printf("\n\n");

        master_sender(A, B, dim, N);
    }

    //    Ablock = matrix_creator(stripSize, N);
    //    Bblock = matrix_creator(stripSize, N);
    //    Cblock = matrix_creator(stripSize, N);
    //
    //    /*receive my part of A and B*/
    //    MPI_Recv(Ablock[0], numElements, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    MPI_Recv(Bblock[0], numElements, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //    printmatrix(stripSize, N, Ablock);

    Ablock = (double *) malloc(N * sizeof (double));
    Bblock = (double *) malloc(N * sizeof (double));
    //Cblock = (double *) malloc(N * sizeof (double));
    Cblock = matrix_creator(dim, dim);

    /*receive my part of A and B*/
    MPI_Recv(Ablock, numElements, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(Bblock, numElements, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /* coords computation */
    coo = coordinate(myrank, numnodes);
    /*printf("My rank is: %d. My coords are: %d, %d\n", myrank, coo[0], coo[1]);*/

    /*communicators creation in order to split blocks that will be used in multiplication*/
    MPI_Comm_split(MPI_COMM_WORLD, coo[0], myrank, &MyComm_row);
    MPI_Comm_split(MPI_COMM_WORLD, coo[1], myrank, &MyComm_col);

    int rsize, csize;
    MPI_Comm_size(MyComm_row, &rsize);
    MPI_Comm_size(MyComm_col, &csize);

    double *rbuf = (double *) malloc(rsize * N * sizeof (double));
    double *cbuf = (double *) malloc(csize * N * sizeof (double));

    //    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Allgather(Ablock, numElements, MPI_DOUBLE, rbuf, numElements, MPI_DOUBLE, MyComm_row);
    MPI_Allgather(Bblock, numElements, MPI_DOUBLE, cbuf, numElements, MPI_DOUBLE, MyComm_col);

    /* restore matrix version */
    /*MAY NOT BE CORRECT*/
    double **AAblock, **BBblock;


    //    printf("\n\n *** Rsize = %d, csize = %d, mult = %d\n", rsize, csize, rsize * N);
    //    printf("NumElements is: %d\n\n", numElements);
    //    printf("Vecotr Ablock:\n");
    //    printvector(N, Ablock);
    //    printf("Vecotr rbuf:\n");
    //    printvector(rsize * N, rbuf);
    //    printf("Vecotr cbuf:\n");
    //    printvector(csize * N, cbuf);


    //    AAblock = devectorizer(rsize, N, rbuf);
    //    BBblock = devectorizer(N, csize, cbuf);

    //test nuovo devettorizzatore
    double ** BBtest = matrix_creator(N, csize);
    x = 0;
    for (i = 0; i < N; i++)
        for (j = 0; j < dim; j++) {
            BBtest[i][j] = cbuf[x];
            x++;
        }
    //    printf("Vecotr cbuf:\n");
    //    printvector(csize * N, cbuf);
    //    printf("BBblock\n");
    //    printmatrix(N, csize, BBblock);
    //    printf("\nBBtest\n");
    //    printmatrix(N, csize, BBtest);
    double ** AAtest = matrix_creator(rsize, N);
    x = 0;
    for (k = 0; k < N; k += dim)
        for (i = 0; i < dim; i++)
            for (j = k; j < k + dim; j++) {
                AAtest[i][j] = rbuf[x];
                x++;
            }

    //    printf("Vecotr rbuf:\n");
    //    printvector(rsize * N, rbuf);
    //    printf("AAblock\n");
    //    printmatrix(rsize, N,AAblock);
    //    printf("\nAAtest\n");
    //    printmatrix(rsize, N, AAtest);

    int l, m;
    //    for (i = 0; i < dim; i++) {
    //        m = 0;
    for (l = 0; l < dim; l++) {
        for (j = 0; j < dim; j++) {
            Cblock[l][j] = 0.0;
            for (k = 0; k < N; k++) {
                Cblock[l][j] += AAtest[l][k] * BBtest[k][j];
            }
            //m++;
        }
    }
    //}

    double * cvector = matrix_vectorizer(dim, dim, Cblock);
    /* master receives from workers  -- note could be done via MPI_Gather */
    double *Carray;

    if (myrank == 0) {
        int gsize;
        MPI_Comm_size(MPI_COMM_WORLD, &gsize);
        Carray = (double *) malloc(N * N * sizeof (double));
    }

    MPI_Gather(cvector, numElements, MPI_DOUBLE, Carray, numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        /* trasform Carray into matrix C */
        /*MAY NOT BE CORRECT*/
        //C = devectorizer(N, N, Carray);

        int r;
        for (r = 0; r < N * N; r++)
            printf("%f ,", Carray[r]);
        printf("\n\n");

        double ** C = matrix_creator(N, N);
        int i, j, x, y, el = 0;
        /*base point (in final matrix) scrolling*/
        /*i.e. 0,0 - 0,2 - 2,0 - 2,2 with N=4 and nproc=4 */
        for (i = 0; i < N; i += dim)
            for (j = 0; j < N; j += dim)
                /*block scorlling, dim: block dimension*/
                for (x = i; x < dim + i; x++)
                    for (y = j; y < dim + j; y++) {
                        printf("x: %d, y: %d, el: %f\n", x,y,Carray[el]);
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
    freematrix(rsize, AAtest);
    freematrix(N, BBtest);
    free(rbuf);
    free(cbuf);
    free(watch);

    MPI_Finalize();
    return 0;
}