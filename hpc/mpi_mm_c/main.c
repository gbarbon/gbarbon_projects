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
    int* coord = (int*) calloc(2, sizeof (int)); //aggiunto (int*)
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    printf("Myrank is %d.Must NOT be 0. Coordinates calculated\n", procNum + 1);
    return coord;
}

int main(int argc, char *argv[]) {
    double **A, **B, **C, *tmp, *tmpA, *tmpB, **Avett, **Bvett;
    double startTime, endTime;
    int numElements, offset, stripSize, myrank, numnodes, N, i, j, k, r, c;

    //commento


    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    MPI_Comm MyComm_row;
    MPI_Comm MyComm_col;
    int* coo;

    N = atoi(argv[1]);
    numnodes = 4;
    //debug
    printf("Printf atoi N: %d\n", N);

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == 0) {
        tmp = (double *) malloc(sizeof (double) * N * N);
        A = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            A[i] = &tmp[i * N];
    } else {
        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
        A = (double **) malloc(sizeof (double *) * N / numnodes);
        for (i = 0; i < N / numnodes; i++)
            A[i] = &tmp[i * N];
    }


    if (myrank == 0) {
        tmp = (double *) malloc(sizeof (double) * N * N);
        B = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            B[i] = &tmp[i * N];
    } else {
        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
        B = (double **) malloc(sizeof (double *) * N / numnodes);
        for (i = 0; i < N / numnodes; i++)
            B[i] = &tmp[i * N];
    }


    if (myrank == 0) {
        tmp = (double *) malloc(sizeof (double) * N * N);
        C = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            C[i] = &tmp[i * N];
    } else {
        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
        C = (double **) malloc(sizeof (double *) * N / numnodes);
        for (i = 0; i < N / numnodes; i++)
            C[i] = &tmp[i * N];
    }

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
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%f ", A[i][j]);
            }
            printf("\n");
        }

        // suddivisione in blocchi della matrice
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
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%f ", Avett[i][j]);
            }
            printf("\n");
        }

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

        for (i = 1; i <= numnodes; i++) {
            MPI_Send(Avett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(Bvett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            printf("rank: %d  ", i);
            for (j = 0; j < N; j++) {
                printf("%f ", Avett[offset][j]);
            }
            printf("offset: %d\n", offset);

            offset += stripSize;

        }
        // si può togliere il for e l'else e usare la scatter
        //MPI_Scatter(Avett, numElements, MPI_DOUBLE, A[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Scatter(Bvett, numElements, MPI_DOUBLE, B[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //debug
        printf("Myrank is %d. Must be 0. Pieces of A and B sent.\n", myrank);
    } else { // receive my part of A and B
        printf("stripSizeR: %d\n", stripSize);
        MPI_Recv(A[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(B[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("rankR: %d  ", myrank);
        for (j = 0; j < stripSize * N; j++) {
            printf("%f ", A[0][j]);
        }
        printf("\n");

        printf("Myrank is %d. Must NOT be 0. Pieces of A and B received.\n", myrank);
        // calcolo delle coordinate
        coo = coordinate(myrank - 1, numnodes);
        printf("Myrank is %d. Must NOT be 0. Coordinates calculated (return to main funct).\n", myrank);
        // creazione communicatori per la condivisione dei blocchi necessari alla moltiplicazione

        printf("Printf coo[0]: %d\n", coo[0]);
        printf("Printf coo[1]: %d\n", coo[1]);

    }

    MPI_Barrier(MPI_COMM_WORLD);

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
    } else {
        for (i = 0; i < N / numnodes; i++) {
            for (j = 0; j < N / numnodes; j++) {
                C[i][j] = 0.0;
            }
        }
    }

    //debug
    printf("Myrank is %d. C initialized\n", myrank);
    double **blocchiA, **blocchiB;
    if (myrank != 0) {
        printf("Myrank is %d. Must NOT be 0. si inizia\n", myrank);


        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
        blocchiA = (double **) malloc(sizeof (double *) * N / numnodes);
        for (i = 0; i < N / numnodes; i++)
            blocchiA[i] = &tmp[i * N];

        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
        blocchiB = (double **) malloc(sizeof (double *) * N / numnodes);
        for (i = 0; i < N / numnodes; i++)
            blocchiB[i] = &tmp[i * N];

        printf("Myrank is %d. Must NOT be 0. AllGather\n", myrank);
        // uso la allGather x ricevere i blocchi e svolgo la moltiplicazione
        MPI_Allgather(A[0], stripSize * N, MPI_DOUBLE, blocchiA[0], stripSize * N, MPI_DOUBLE, MyComm_row);
        MPI_Allgather(B[0], stripSize * N, MPI_DOUBLE, blocchiB[0], stripSize * N, MPI_DOUBLE, MyComm_col);


        printf("rankM: %d  ", myrank);
        for (j = 0; j < stripSize * N; j++) {
            printf("%f ", blocchiA[0][j]);
        }
        printf("\n");
        printf("Myrank is %d. Must NOT be 0. Multiplies\n", myrank);
        // do the work
        /*for (i = 0; i <= stripSize; i++) {
            for (j = 0; j <= stripSize; j++) {
                for (k = 0; k <= stripSize; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                    C[i][j] += blocchiA[i][k] * blocchiB[k][j];
                }
            }
        }*/
        //debug
        printf("Myrank is %d. Must NOT be 0. Work done!!\n", myrank);
    }

    // master receives from workers  -- note could be done via MPI_Gather
    /*if (myrank == 0) {
        printf("Myrank is %d. Must be 0. Calculate\n", myrank);
        offset = stripSize;
        numElements = stripSize * N;
        for (i = 1; i <= numnodes; i++) {
            MPI_Recv(C[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            offset += stripSize;
        }
        //debug
        printf("Myrank is %d. Must be 0. Pieces received from workers\n", myrank);
    } else { // send my contribution to C
        MPI_Send(C[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
        //debug

        printf("Myrank is %d. Must NOT be 0, I am a worker. My contribution has been sent\n", myrank);

    }*/

    // stop timer
    if (myrank == 0) {
        endTime = MPI_Wtime();
        printf("Time is %f\n", endTime - startTime);
        free(A);
        free(B);
        free(C);
        free(tmp);
        free(tmpA);
        free(tmpB);
        free(Avett);
        free(Bvett);
    }


    // print out matrix here, if I'm the master
    if (myrank == 0 && N < 10) {
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                printf("%f ", C[i][j]);
            }
            printf("\n");
        }
    }


    /*free(blocchiA);
    free(blocchiB);
    free(coord);*/

    MPI_Finalize();
    return 0;
}



