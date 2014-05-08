/* 
 * File:   mult_pf.c
 * Author: asus
 *
 * Created on 24 marzo 2014, 10.53
 */

#include "header.h"
#include "MpiStopwatch.h"
#include "inout.h"
#define TAG 13

int main(int argc, char** argv) {

    double **A, **B, **C, *rigaA, *ris, *Bvett;
    int recv, indexR, index = 0, myrank, numnodes, N, i, j;

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    N = atoi(argv[1]);

    /*I\O*/
    int inout_bool = atoi(argv[2]);
    char * homePath = getenv("HOME"); /*homepath*/

    // heavy load f(A) abilitation
    int load_bool = atoi(argv[3]);

    /*CSV file support*/
    char *op = "farm", final[256];
    int myid;

    /*stopwatch*/
    Stopwatch watch = StopwatchCreate();

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == 0) {
        printf("Printf atoi N: %d\n", N);
        printf("Printf numnodes: %d\n", numnodes);

        /*start timer*/
        StopwatchStart(watch);

        /*input type evaluation*/
        if (inout_bool == 0) { /*no input, so randomly generated matrix*/

            /* matrix creation */
            A = matrix_creator(N, N);
            B = matrix_creator(N, N);

            /*init matrices with random values*/
            simple_matrix_init(A, N);
            simple_matrix_init(B, N);

        } else if (inout_bool == 1) {
            /*input filename generation*/
            char infile[256];
            snprintf(infile, sizeof infile, "%s/hpc_temp/hpc_input/mat%d.csv", homePath, N);

            /* matrix loading */
            A = matrix_loader(infile);
            B = matrix_loader(infile);
        } else {
            printf("Error on input output value!!!\n");
            return 0;
        }

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
                    if (load_bool == 0)
                        ris[j] += rigaA[i] * B[i][j];
                    else
                        ris[j] += heavy(rigaA[i], load_bool) * B[i][j];
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

        if (recv == N) {

            /*output type evaluation*/
            if (inout_bool == 0) { /*no output, so no result (or print)*/
                /*printmatrix(N, N, res);*/
            } else {
                char outfile[256];
                snprintf(outfile, sizeof outfile, "%s/hpc_temp/hpc_output/%s_dim%d_nproc%d.csv", homePath, op, N, numnodes - 1);
                matrix_writer(N, C, outfile);
            }

            /*stopwatch stop*/
            StopwatchStop(watch);
            StopwatchPrintWithComment("Master total time: %f\n", watch);
            myid = (int) MPI_Wtime(); /*my id generation*/
            snprintf(final, sizeof final, "%d,%s%s,%d,%d,%d,%d", myid, op, OPTI, numnodes - 1, N, inout_bool, load_bool); /*final string generation*/
            StopwatchPrintToFile(final, watch);

            freematrix(N, A);
            freematrix(N, C);
        }
    }

    freematrix(N, B);
    free(Bvett);
    free(watch);

    if (myrank != 0) {
        free(rigaA);
        free(ris);
    }


    MPI_Finalize();

    return 0;
}





