/* 
 * File:   mpi_mm.cpp
 * Author: jian
 *
 * Created on 11 marzo 2014, 17.22
 */


// MPI matrix matrix multiplication

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TAG 13

int* coordinate(int procNum, int totalProc) {
    int* coord = calloc(2, sizeof (int));
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    return coord;
}

int main(int argc, char *argv[]) {
    double **A, **B, **C, *tmp, *tmpA, *tmpB, **Avett, **Bvett;
    double startTime, endTime;
    int numElements, offset, stripSize, myrank, numnodes, N, i, j, k, r, c;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    N = atoi(argv[1]);

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

    if (myrank == 0) {
        // initialize A and B
        double w = 0.0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                A[i][j] = k;
                B[i][j] = k;
                w = w + 1.0;
            }
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

    }

    // start timer
    if (myrank == 0) {
        startTime = MPI_Wtime();
    }

    stripSize = N / numnodes;

    // send each node its piece of A and B
    if (myrank == 0) {
        offset = stripSize;
        numElements = stripSize * N;

        for (i = 1; i < numnodes; i++) {
            MPI_Send(Avett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(Bvett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            offset += stripSize;
        }
        // si può togliere il for e l'else e usare la scatter
        //MPI_Scatter(Avett, numElements, MPI_DOUBLE, A[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        //MPI_Scatter(Bvett, numElements, MPI_DOUBLE, B[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    } else { // receive my part of A and B
        MPI_Recv(A[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(B[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // calcolo delle coordinate
        int* coo = coordinate(myrank - 1, numnodes);

        // creazione communicatori per la condivisione dei blocchi necessari alla moltiplicazione
        int MPI_Comm_split(MPI_COMM_WORLD, coo[0], myrank, &MyComm_row);
        int MPI_Comm_split(MPI_COMM_WORLD, coo[1], myrank, &MyComm_col);
    }

    // Let each process initialize C to zero 
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            C[i][j] = 0.0;
        }
    }

    if (myrank != 0) {

        int **blocchiA, **blocchiB;

        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
        blocchiA = (double **) malloc(sizeof (double *) * N / numnodes);
        for (i = 0; i < N / numnodes; i++)
            blocchiA[i] = &tmp[i * N];

        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
        blocchiB = (double **) malloc(sizeof (double *) * N / numnodes);
        for (i = 0; i < N / numnodes; i++)
            blocchiB[i] = &tmp[i * N];

        // uso la allGather x ricevere i blocchi e svolgo la moltiplicazione
        MPI_Allgather(A[0], stripSize * N, MPI_DOUBLE, blocchiA, stripSize * N, MPI_DOUBLE, MyComm_row);
        MPI_Allgather(B[0], stripSize * N, MPI_DOUBLE, blocchiB, stripSize * N, MPI_DOUBLE, MyComm_col);

        // do the work
        for (i = 0; i <= stripSize; i++) {
            for (j = 0; j <= stripSize; j++) {
                for (k = 0; k <= stripSize; k++) {
                    C[i][j] += A[i][k] * B[k][j];
                    C[i][j] += blocchiA[i][k] * blocchiB[k][j];
                }
            }
        }


    }

    // master receives from workers  -- note could be done via MPI_Gather
    if (myrank == 0) {
        offset = stripSize;
        numElements = stripSize * N;
        for (i = 1; i < numnodes; i++) {
            MPI_Recv(C[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            offset += stripSize;
        }
    } else { // send my contribution to C
        MPI_Send(C[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
    }

    // stop timer
    if (myrank == 0) {
        endTime = MPI_Wtime();
        printf("Time is %f\n", endTime - startTime);
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

    MPI_Finalize();
    return 0;
}


