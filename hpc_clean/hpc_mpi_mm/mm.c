
/* 
 * File:   mpi_mm.c
 * Author: 
 *
 * Created on 12 marzo 2014, 15.20
 */

/* MPI matrix matrix multiplication */

#include "header.h"

int main(int argc, char *argv[]) {
    //    double **A, **B, **C, *tmp, *tmpA, *tmpB, *tmpC, **Avett, **Bvett;
    //    double **Ablock, **Bblock, **Cblock;
    //    double startTime, endTime;
    //    int numElements, offset, stripSize, N, i, j, k, r, c;

    /*MPI variables*/
    MPI_Comm MyComm_row, MyComm_col;
    MPI_Status status;
    int myrank, numnodes;

    /*Test variables*/
    int ind_split, req_tag = 0, ans_tag = 2;
    char message[100];


    /*MPI initialization*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    /*show who I am*/
    printf("I'm process %d\n", myrank);

    if (myrank == 0) {
        /*split matrix in pieces*/
        /*send matrix pieces*/
        for (ind_split = 1; ind_split <= numnodes - 1; ind_split++) {
            printf("Debug: dentro il ciclo %d\n", ind_split);
            sprintf(message, "This message for processor number %d\n", ind_split);
            MPI_Send(message, strlen(message) + 1, MPI_CHAR, ind_split, req_tag, MPI_COMM_WORLD);
        }
    }

        /*barrier to wait that all nodes receive the message*/
        /*MPI_Barrier(MPI_COMM_WORLD);*/
        /*substitued with a synchronous send*/
    else {
        /*the process receive it's part*/
        MPI_Recv(message, 100, MPI_CHAR, 0, req_tag, MPI_COMM_WORLD, &status);
        /*work*/
        sprintf(message, "%s COMPUTED BY PROCESSOR NUMBER: %d\n", message, myrank);
        /*send work back to master*/
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, ans_tag, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (myrank == 0) {
        //MPI_Reduce();
        printf("*** FINAL PRINT ***\n");
        for (ind_split = 1; ind_split <= numnodes - 1; ind_split++) {
            MPI_Recv(message, 100, MPI_CHAR, ind_split, ans_tag, MPI_COMM_WORLD, &status);
            printf("Printing message: %s\n", message);
        }
    }

    //    int* coo;
    //
    /*variables for number of processors*/
    /*unused if we use parameter -n on mpiexec, program will take n. of processor and handle with 'numnodes' var*/
    //N = atoi(argv[1]);
    //numnodes = 4;
    //
    //    /*debug*/
    //    printf("Myrank is %d.\n", myrank);
    //
    //    /* allocate A, B, and C --- note that you want these to be contiguously allocated.  Workers need less memory allocated. */
    //
    //    if (myrank == 0) {
    //
    //        /*A, B, C allocation*/
    //        /*tmp = (double *) malloc(sizeof (double) * N * N);
    //        A = (double **) malloc(sizeof (double *) * N);
    //        for (i = 0; i < N; i++)
    //            A[i] = &tmp[i * N];
    //
    //        tmp = (double *) malloc(sizeof (double) * N * N);
    //        B = (double **) malloc(sizeof (double *) * N);
    //        for (i = 0; i < N; i++)
    //            B[i] = &tmp[i * N];
    //
    //        tmp = (double *) malloc(sizeof (double) * N * N);
    //        C = (double **) malloc(sizeof (double *) * N);
    //        for (i = 0; i < N; i++)
    //            C[i] = &tmp[i * N];*/
    //
    //        //debug
    //        printf("Myrank is %d. Should be 0. A,B,C allocated\n", myrank);
    //
    //        // initialize A and B
    //        double w = 0.0;
    //        for (i = 0; i < N; i++) {
    //            for (j = 0; j < N; j++) {
    //                A[i][j] = w;
    //                B[i][j] = w;
    //                w = w + 1.0;
    //            }
    //        }
    //
    //        //debug
    //        printf("Myrank is %d. Should be 0. A,B initialized\n", myrank);
    //
    //        // suddivisione in blocchi della matrice
    //        tmpA = (double *) malloc(sizeof (double) * N * N);
    //        tmpB = (double *) malloc(sizeof (double) * N * N);
    //        Avett = (double **) malloc(sizeof (double *) * N);
    //        Bvett = (double **) malloc(sizeof (double *) * N);
    //
    //        int k = 0;
    //        for (i = 0; i < N; i = i + N / (numnodes / 2)) {
    //            for (j = 0; j < N; j = j + N / (numnodes / 2)) {
    //                for (r = i; r < i + N / (numnodes / 2); r++) {
    //                    for (c = j; c < j + N / (numnodes / 2); c++) {
    //                        tmpA[k] = A[r][c];
    //                        tmpB[k] = B[r][c];
    //                        k++;
    //                    }
    //                }
    //            }
    //        }
    //
    //
    //        for (i = 0; i < N; i++) {
    //            Avett[i] = &tmpA[i * N];
    //            Bvett[i] = &tmpB[i * N];
    //        }
    //
    //        //debug
    //        printf("Myrank is %d. Should be 0. A,B divided\n", myrank);
    //
    //    } else {
    //        /*tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        A = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            A[i] = &tmp[i * N];*/
    //
    //        tmpA = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        Ablock = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            Ablock[i] = &tmpA[i * N];
    //
    //        /*tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        B = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            B[i] = &tmp[i * N];*/
    //
    //        tmpB = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        Bblock = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            Bblock[i] = &tmpB[i * N];
    //
    //        /*tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        C = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            C[i] = &tmp[i * N]; */
    //        tmpC = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        Cblock = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            Cblock[i] = &tmpC[i * N];
    //
    //        //debug
    //        printf("Myrank is %d. Should NOT be 0. A,B,C allocated\n", myrank);
    //    }
    //
    //    // start timer
    //    if (myrank == 0) {
    //        startTime = MPI_Wtime();
    //    }
    //
    //    stripSize = N / numnodes;
    //
    //    // send each node its piece of A and B
    //    if (myrank == 0) {
    //        offset = stripSize;
    //        numElements = stripSize * N;
    //
    //        for (i = 1; i < numnodes; i++) {
    //            MPI_Send(Avett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
    //            MPI_Send(Bvett[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
    //            offset += stripSize;
    //        }
    //
    //        //debug
    //        printf("Myrank is %d. Should be 0. Pieces sent\n", myrank);
    //
    //        // si può togliere il for e l'else e usare la scatter
    //        //MPI_Scatter(Avett, numElements, MPI_DOUBLE, A[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //        //MPI_Scatter(Bvett, numElements, MPI_DOUBLE, B[0], numElements, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    //    } else { // receive my part of A and B
    //        /*
    //        MPI_Recv(A[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //        MPI_Recv(B[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //        */
    //        MPI_Recv(Ablock[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //        MPI_Recv(Bblock[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //
    //        //debug
    //        printf("Myrank is %d. Should NOT be 0. Pieces received\n", myrank);
    //
    //        // calcolo delle coordinate
    //        int* coo = coordinate(myrank, numnodes);
    //
    //        //debug
    //        printf("Myrank is %d. Should NOT be 0. Coords computed\n", myrank);
    //        printf("Printf coo[0]: %d\n", coo[0]);
    //        printf("Printf coo[1]: %d\n", coo[1]);
    //
    //        // creazione communicatori per la condivisione dei blocchi necessari alla moltiplicazione
    //
    //        MPI_Comm_split(MPI_COMM_WORLD, coo[0], myrank, &MyComm_row);
    //        MPI_Comm_split(MPI_COMM_WORLD, coo[1], myrank, &MyComm_col);
    //
    //        //debug
    //        printf("Myrank is %d. Should NOT be 0. Communicators created\n", myrank);
    //    }
    //
    //    // Let each process initialize C to zero 
    //    for (i = 0; i < N; i++) {
    //        for (j = 0; j < N; j++) {
    //            C[i][j] = 0.0;
    //        }
    //    }
    //
    //    //debug
    //    printf("Myrank is %d. C initialized\n", myrank);
    //
    //    if (myrank != 0) {
    //
    //        double **blocchiA, **blocchiB;
    //
    //        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        blocchiA = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            blocchiA[i] = &tmp[i * N];
    //
    //        tmp = (double *) malloc(sizeof (double) * N * N / numnodes);
    //        blocchiB = (double **) malloc(sizeof (double *) * N / numnodes);
    //        for (i = 0; i < N / numnodes; i++)
    //            blocchiB[i] = &tmp[i * N];
    //
    //        // uso la allGather x ricevere i blocchi e svolgo la moltiplicazione
    //        MPI_Allgather(A[0], stripSize * N, MPI_DOUBLE, blocchiA, stripSize * N, MPI_DOUBLE, MyComm_row);
    //        MPI_Allgather(B[0], stripSize * N, MPI_DOUBLE, blocchiB, stripSize * N, MPI_DOUBLE, MyComm_col);
    //
    //        // do the work
    //        for (i = 0; i <= stripSize; i++) {
    //            for (j = 0; j <= stripSize; j++) {
    //                for (k = 0; k <= stripSize; k++) {
    //                    C[i][j] += A[i][k] * B[k][j];
    //                    C[i][j] += blocchiA[i][k] * blocchiB[k][j];
    //                }
    //            }
    //        }
    //
    //
    //    }
    //
    //    // master receives from workers  -- note could be done via MPI_Gather
    //    if (myrank == 0) {
    //        offset = stripSize;
    //        numElements = stripSize * N;
    //        for (i = 1; i < numnodes; i++) {
    //            MPI_Recv(C[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //            offset += stripSize;
    //        }
    //
    //        // stop timer
    //
    //        endTime = MPI_Wtime();
    //        printf("Time is %f\n", endTime - startTime);
    //
    //    } else { // send my contribution to C
    //        MPI_Send(C[0], stripSize * N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
    //    }
    //
    //    // stop timer
    //    if (myrank == 0) {
    //        //endTime = MPI_Wtime();
    //        //printf("Time is %f\n", endTime - startTime);
    //        free(A);
    //        free(B);
    //        free(C);
    //        free(tmp);
    //        free(tmpA);
    //        free(tmpB);
    //        free(tmpC);
    //        free(Avett);
    //        free(Bvett);
    //        free(Ablock);
    //        free(Bblock);
    //        free(Cblock);
    //    }
    //
    //    /*print out matrix here, if I'm the master*/
    //    /*if (myrank == 0 && N < 10) {
    //        printmatrix(N, C);
    //    }*/
    //
    //    //free all 
    MPI_Finalize();
    return 0;
}



