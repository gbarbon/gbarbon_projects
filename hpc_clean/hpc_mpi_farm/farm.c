/* 
 * File:   mult_pf.c
 * Author: asus
 *
 * Created on 24 marzo 2014, 10.53
 */

#include "header.h"
#define TAG 13

int main(int argc, char** argv) {

    double **A, **B, **C, *tmp, *rigaA, *ris, *Bvett;
    double startTime, endTime;
    int recv, indexR, index, myrank, numnodes, N, i, j, k, r, c;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    N = atoi(argv[1]);

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == 0) {
        printf("Printf atoi N: %d\n", N);
        printf("Printf numnodes: %d\n", numnodes);

        A = matrix_creator(N, N);

        B = matrix_creator(N, N);

        C = matrix_creator(N, N);

        Bvett = (double *) malloc(sizeof (double) * N * N);

        printf("Myrank is %d. A,B,C allocated\n", myrank);
    }

    //debug

    if (myrank != 0) {
        Bvett = (double *) malloc(sizeof (double) * N * N);
        rigaA = (double *) malloc(N * sizeof (double));
        ris = (double *) malloc(N * sizeof (double));
    }

    if (myrank == 0) {
        // initialize A and B
        simple_matrix_init(A, N);
        simple_matrix_init(B, N);

        Bvett = matrix_vectorizer(N, N, B);

        for (i = 0; i < numnodes - 1; i++) {
            MPI_Send(Bvett, N*N, MPI_DOUBLE, i + 1, TAG, MPI_COMM_WORLD);
        }
        printf("Myrank is %d. B inviated\n", myrank);
    } else {
        MPI_Recv(Bvett, N*N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        B = devectorizer(N, N, Bvett);
        printf("Myrank is %d. B received\n", myrank);
    }

    if (myrank == 0) {

        // start timer
        startTime = MPI_Wtime();

        // invio di una riga di A ad ogni nodo slave
        for (i = 0; i < numnodes - 1; i++) {
            if (index < N && index != -1) {
                MPI_Send(&index, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
                MPI_Send(A[index], N, MPI_DOUBLE, i + 1, i + 1, MPI_COMM_WORLD);
                index = index + 1;
            }// se ci sono più nodi che righe mando -1
            else {
                index = -1;
                MPI_Send(&index, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
            }
        }
        printf("Myrank is %d. Send effettuate\n", myrank);

    }

    if (myrank != 0) {
        // ricevo l'indice e se è != -1 ricevo la riga di A ed eseguo la moltiplicazione
        MPI_Recv(&index, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        while (index != -1) {
            MPI_Recv(rigaA, N, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            for (j = 0; j < N; j++) ris[j] = 0.0;

            for (j = 0; j < N; j++) {
                for (i = 0; i < N; i++) {
                    ris[j] += rigaA[i] * B[i][j];
                }
            }
            printf("Myrank is %d. Recv effettuate\n", myrank);

            // effettuata la moltiplicazione invio al master il risultato
            MPI_Send(&index, 1, MPI_INT, 0, TAG, MPI_COMM_WORLD);
            MPI_Send(ris, N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD);
            printf("Myrank is %d. Mult effettuate\n", myrank);

            MPI_Recv(&index, 1, MPI_INT, 0, myrank, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }


    }

    if (myrank == 0) {
        recv = 0;
        // ricevo dagli slave le righe calcolate
        while (recv < N) {
            for (i = 0; i < numnodes - 1; i++) {
                if (recv < N) {
                    MPI_Recv(&indexR, 1, MPI_INT, i + 1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(C[indexR], N, MPI_DOUBLE, i + 1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    recv++;
                }
                // se ci sono ancora righe da inviare le invio al processo appena liberatosi
                if (index < N && index != -1) {
                    MPI_Send(&index, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
                    MPI_Send(A[index], N, MPI_DOUBLE, i + 1, i + 1, MPI_COMM_WORLD);
                    index = index + 1;
                }// se ci sono più nodi che righe mando -1
                else {
                    index = -1;
                    MPI_Send(&index, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
                }
            }
        }

        // stop timer
        endTime = MPI_Wtime();
        printf("Time is %f\n", endTime - startTime);


        if (recv == N) {
            printmatrix(N, N, C);

            freematrix(N, A);
            freematrix(N, C);
        }
    }

    freematrix(N, B);
    free(Bvett);

    if (myrank != 0) {
        free(rigaA);
        free(ris);
    }


    MPI_Finalize();

    return 0;
}





