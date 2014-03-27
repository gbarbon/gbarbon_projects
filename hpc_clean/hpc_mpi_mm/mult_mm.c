/*
 * File:   mpi_mm.c
 * Author: jian
 *
 * Created on 12 marzo 2014, 15.20
 */

// MPI matrix matrix multiplication

#include <mpi.h>


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAG 13

int* coordinate(int procNum, int totalProc) {
    int* coord = (int*) calloc(2, sizeof (int));
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    printf("Myrank is %d.Must NOT be 0. Coordinates calculated\n", procNum);
    return coord;
}

int main(int argc, char *argv[]) {
    double **A, **B, **C, *tmp, *tmpA, *tmpB, *tmpC, **Avett, **Bvett;
    double **Ablock, **Bblock, **Cblock;
    double startTime, endTime;
    int numElements, offset, stripSize, myrank, numnodes, N, i, j, k, r, c;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    MPI_Comm MyComm_row;
    MPI_Comm MyComm_col;
    int* coo;

    N = atoi(argv[1]);
    //numnodes = 4;
    //debug
    printf("Printf atoi N: %d\n", N);
    printf("Printf numnodes: %d\n", numnodes);

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == 0) {
        tmp = (double *) malloc(sizeof (double) * N * N);
        A = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            A[i] = &tmp[i * N];
    }
    tmpA = (double *) malloc(sizeof (double) * N * N / numnodes);
    Ablock = (double **) malloc(sizeof (double *) * N / numnodes);
    for (i = 0; i < N / numnodes; i++)
        Ablock[i] = &tmpA[i * N];



    if (myrank == 0) {
        tmp = (double *) malloc(sizeof (double) * N * N);
        B = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            B[i] = &tmp[i * N];
    }

    tmpB = (double *) malloc(sizeof (double) * N * N / numnodes);
    Bblock = (double **) malloc(sizeof (double *) * N / numnodes);
    for (i = 0; i < N / numnodes; i++)
        Bblock[i] = &tmpB[i * N];



    if (myrank == 0) {
        tmp = (double *) malloc(sizeof (double) * N * N);
        C = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            C[i] = &tmp[i * N];
    }
    tmpC = (double *) malloc(sizeof (double) * N * N / numnodes);
    Cblock = (double **) malloc(sizeof (double *) * N / numnodes);
    for (i = 0; i < N / numnodes; i++)
        Cblock[i] = &tmpC[i * N];


    //debug
    printf("Myrank is %d. A,B,C allocated\n", myrank);

    if (myrank == 0) {
        // initialize A and B
        double w = 0.0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                A[i][j] = w;
                B[i][j] = w;
                w = w + 1.0;
            }
        }
        /*for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%f ", A[i][j]);
            }
            printf("\n");
        }*/

        // suddivisione in blocchi delle matrici
        tmpA = (double *) malloc(sizeof (double) * N * N);
        tmpB = (double *) malloc(sizeof (double) * N * N);
        Avett = (double **) malloc(sizeof (double *) * N);
        Bvett = (double **) malloc(sizeof (double *) * N);

        // suddivisione in blocchi della matrice A
        int k = 0;
        for (i = 0; i < N; i = i + N / (numnodes / 2)) {
            for (j = 0; j < N; j = j + N / (numnodes / 2)) {
                for (r = i; r < i + N / (numnodes / 2); r++) {
                    for (c = j; c < j + N / (numnodes / 2); c++) {
                        tmpA[k] = A[r][c];
                        k++;
                    }
                }
            }
        }

        // suddivisione in blocchi della matrice B
        k = 0;
        for (i = 0; i < N; i = i + N / (numnodes / 2)) {
            for (j = 0; j < N; j = j + N / (numnodes / 2)) {
                for (r = j; r < j + N / (numnodes / 2); r++) {
                    for (c = i; c < i + N / (numnodes / 2); c++) {
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
        /*for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%f ", Avett[i][j]);
            }
            printf("\n");
        }*/

    }

    // start timer
    if (myrank == 0) {
        startTime = MPI_Wtime();
    }

    stripSize = N / numnodes;

    // send each node its piece of A and B
    if (myrank == 0) {

        offset = 0;
        numElements = stripSize * N;
        printf("stripSizeS: %d\n", stripSize);

        for (i = 0; i < numnodes; i++) {
            MPI_Send(Avett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(Bvett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);

            /*printf("rank: %d  ", i);
            for (j = 0; j < N; j++) {
                printf("%f ", Avett[offset][j]);
            }
            printf("offset: %d\n", offset);*/

            offset += stripSize;

        }
        // si può togliere il for e l'else e usare la scatter
        //MPI_Scatter(Avett, numElements, MPI_DOUBLE, A[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Scatter(Bvett, numElements, MPI_DOUBLE, B[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //debug
        printf("Myrank is %d. Must be 0. Pieces of A and B sent.\n", myrank);
    }
    // receive my part of A and B
    printf("stripSizeR: %d\n", stripSize);
    MPI_Recv(Ablock[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(Bblock[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    /*printf("rankR: %d  ", myrank);
    for (j = 0; j < stripSize * N; j++) {
        printf("%f ", Ablock[0][j]);
    }
    printf("\n");*/

    printf("Myrank is %d. Must NOT be 0. Pieces of A and B received.\n", myrank);
    // calcolo delle coordinate
    coo = coordinate(myrank, numnodes);
    printf("Myrank is %d. Must NOT be 0. Coordinates calculated (return to main funct).\n", myrank);

    // creazione communicatori per la condivisione dei blocchi necessari alla moltiplicazione

    printf("Printf coo[0]: %d\n", coo[0]);
    printf("Printf coo[1]: %d\n", coo[1]);

    MPI_Comm_split(MPI_COMM_WORLD, coo[0], myrank, &MyComm_row);
    MPI_Comm_split(MPI_COMM_WORLD, coo[1], myrank, &MyComm_col);

    //debug
    printf("Myrank is %d. Must NOT be 0. Communicators created.\n", myrank);



    // Let each process initialize C to zero
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


    //debug
    printf("Myrank is %d. C initialized\n", myrank);

    int rsize, csize;
    double *rbuf;
    double *cbuf;


    printf("Myrank is %d. Must NOT be 0. AllGatherInizio\n", myrank);

    MPI_Comm_size(MyComm_row, &rsize);
    MPI_Comm_size(MyComm_col, &csize);

    printf("Myrank is %d. Must NOT be 0. rsize: %d csize: %d\n", myrank, rsize, csize);

    rbuf = (double *) malloc(rsize * N * sizeof (double));
    cbuf = (double*) malloc(csize * N * sizeof (double));

    printf("Myrank is %d. Must NOT be 0. AllGather\n", myrank);

    MPI_Allgather(Ablock[0], stripSize * N, MPI_DOUBLE, rbuf, stripSize * N, MPI_DOUBLE, MyComm_row);
    MPI_Allgather(Bblock[0], stripSize * N, MPI_DOUBLE, cbuf, stripSize * N, MPI_DOUBLE, MyComm_col);

    printf("Myrank is %d. Must NOT be 0. AllGatherFatta\n", myrank);

    /*printf("rankM: %d  ", myrank);
    for (j = 0; j < rsize * stripSize * N; j++) {
        printf("%f ", rbuf[j]);
    }
    printf("\n");

    printf("rankM: %d  ", myrank);
    for (j = 0; j < csize * stripSize * N; j++) {
        printf("%f ", cbuf[j]);
    }
    printf("\n");*/

    // ripristina la versione a matrice
    double *tmpAA, *tmpBB;
    double **AAblock, **BBblock;

    tmpAA = (double *) malloc(rsize * N * sizeof (double));
    AAblock = (double **) malloc(sizeof (double *) * rsize);
    for (i = 0; i < rsize; i++)
        AAblock[i] = &tmpAA[i * N];

    tmpBB = (double *) malloc(csize * N * sizeof (double));
    BBblock = (double **) malloc(sizeof (double *) * csize);
    for (i = 0; i < csize; i++)
        BBblock[i] = &tmpBB[i * N];


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

    /*printf("rankSA: %d  \n", myrank);
    for (i = 0; i < N / (numnodes / 2); i++) {
        for (j = 0; j < N; j++) {
            printf("%f ", AAblock[i][j]);
        }
        printf("\n");
    }
    printf("rankSB: %d  \n", myrank);
    for (i = 0; i < N / (numnodes / 2); i++) {
        for (j = 0; j < N; j++) {
            printf("%f ", BBblock[i][j]);
        }
        printf("\n");
    }

    printf("\n");*/
    printf("Myrank is %d. Must NOT be 0. Multiplies\n", myrank);
    // do the work

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
    //debug
    printf("Myrank is %d. Must NOT be 0. Work done!!\n", myrank);

    /*printf("rankCB: %d  \n", myrank);
    for (i = 0; i < N / numnodes; i++) {
        for (j = 0; j < stripSize * N; j++) {
            printf("%f ", Cblock[i][j]);
        }
        printf("\n");
    }

    printf("\n");*/


    // master receives from workers  -- note could be done via MPI_Gather
    int gsize;
    double *Carray;

    if (myrank == 0) {
        MPI_Comm_size(MPI_COMM_WORLD, &gsize);
        Carray = (double *) malloc(gsize * N * sizeof (double));
    }

    printf("Gather C\n");
    MPI_Gather(Cblock[0], stripSize * N, MPI_DOUBLE, Carray, stripSize * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    printf("C pieces received from workers\n");

    // stop timer
    if (myrank == 0) {

        // trasformo l'array Carray nella Matrice C
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

        endTime = MPI_Wtime();
        printf("Time is %f\n", endTime - startTime);

        // print out matrix here, if I'm the master
        if (N < 10) {
            printf("Stampa C\n");
            for (i = 0; i < N; i++) {
                for (j = 0; j < N; j++) {
                    printf("%f ", C[i][j]);
                }
                printf("\n");
            }
        }

        free(A);
        free(B);
        free(C);
        free(tmp);
        free(Avett);
        free(Bvett);
        free(Carray);
    }
    
    free(tmpA);
    free(tmpB);
    free(tmpC);
    free(Ablock);
    free(Bblock);
    free(Cblock);
    free(coo);
    free(rbuf);
    free(cbuf);
    free(tmpAA);
    free(tmpBB);
    free(AAblock);
    free(BBblock);

    MPI_Finalize();
    return 0;
}



