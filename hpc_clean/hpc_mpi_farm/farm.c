/* 
 * File:   mult_pf.c
 * Author: asus
 *
 * Created on 24 marzo 2014, 10.53
 */

#include <mpi.h>


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define TAG 13

/*
 * 
 */

double* mat2array(double **Mat, int n, double* vett) {
    int i, j, k;
    k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            vett[k] = Mat[i][j];
            k++;
        }
    }
    return vett;
}

double** array2mat(double* vett, int n, double **Mat) {
    int i, j, k;
    k = 0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            Mat[i][j] = vett[k];
            k++;
        }
    }
    return Mat;
}

void stampaMat(double **Mat, int n) {
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            printf("%f ", Mat[i][j]);
        }
        printf("\n");
    }
}

int main(int argc, char** argv) {

    double **A, **B, **C, *tmp, *rigaA, *ris, *Bvett;
    double startTime, endTime;
    int recv, indexR, index, myrank, numnodes, N, i, j, k, r, c;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    N = atoi(argv[1]);
    //numnodes = 4;
    //debug


    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == 0) {
        printf("Printf atoi N: %d\n", N);
        printf("Printf numnodes: %d\n", numnodes);
        tmp = (double *) malloc(sizeof (double) * N * N);
        A = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            A[i] = &tmp[i * N];

        tmp = (double *) malloc(sizeof (double) * N * N);
        B = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            B[i] = &tmp[i * N];

        tmp = (double *) malloc(sizeof (double) * N * N);
        C = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            C[i] = &tmp[i * N];

        Bvett = (double *) malloc(sizeof (double) * N * N);

        printf("Myrank is %d. A,B,C allocated\n", myrank);
    }

    //debug

    if (myrank != 0) {
        tmp = (double *) malloc(sizeof (double) * N * N);
        B = (double **) malloc(sizeof (double *) * N);
        for (i = 0; i < N; i++)
            B[i] = &tmp[i * N];

        Bvett = (double *) malloc(sizeof (double) * N * N);

        rigaA = (double *) malloc(N * sizeof (double));
        ris = (double *) malloc(N * sizeof (double));

        printf("Myrank is %d. B allocated\n", myrank);
    }

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

        Bvett = mat2array(B, N, Bvett);

        for (i = 0; i < numnodes - 1; i++) {
            MPI_Send(Bvett, N*N, MPI_DOUBLE, i + 1, TAG, MPI_COMM_WORLD);
        }
        printf("Myrank is %d. B inviated\n", myrank);
    } else {
        MPI_Recv(Bvett, N*N, MPI_DOUBLE, 0, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        B = array2mat(Bvett, N, B);
        printf("Myrank is %d. B received\n", myrank);
    }

    if (myrank == 0) {
        // invio di una riga di A ad ogni nodo slave
        for (i = 0; i < numnodes - 1; i++) {
            if (index < N && index != -1) {
                MPI_Send(&index, 1, MPI_INT, i + 1, i + 1, MPI_COMM_WORLD);
                //rigaA = getRiga(A, *indice);
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

    //while (recv < N) {
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
                    printf("indexR= %d\n", indexR);
                    //MPI_Recv(C[indexR], N, MPI_DOUBLE, i + 1, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
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
                printf("Myrank is %d. recv= %d\n", myrank, recv);
            }
        }

        if (recv == N) {
        printf("Stampa C\n");
        stampaMat(C, N);

        free(A);
        free(C);
        }
    }
    //}

    free(B);
    free(tmp);
    free(Bvett);

    if (myrank != 0) {
        free(rigaA);
        free(ris);
    }

    printf("Myrank is %d. Fine while\n", myrank);
    MPI_Finalize();
    printf("Myrank is %d. Finalize\n", myrank);
    return 0;
}





