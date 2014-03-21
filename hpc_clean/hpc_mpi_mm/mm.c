
/* 
 * File:   mm.c
 * Author: 
 *
 * Created on 12 marzo 2014, 15.20
 */

/* MPI matrix matrix multiplication */

#include "header.h"

/* 
 * Compute multiplication
 */
double** mult(double** rows, double** cols, int n, int offset) {
    int i, j, k;
    double** res; /*data structure for resulting matrix*/

    res = matrix_creator(offset, offset);
    for (i = 0; i <= n; i++) {
        for (j = 0; j <= n; j++) {
            for (k = 0; k <= offset; k++) {
                res[i][j] += rows[i][k] * cols[k][j];
            }
        }
    }

    return res;
}

/*
 * Description: send rows and cols of A and B to the interested worker
 */
int master_sender(double** A, double** B, int offset, int n) {
    int i, j, worker = 0;
    for (j = 0; j < n; j += offset)
        for (i = 0; i < n; i += offset) {
            worker++;
            printf("master_sender for worker: %d, offset: %d, n: %d \n\n", worker, offset, n);
            printmatrix(n,n, A);
            printf("\n\n");
            printmatrix(offset,n, &A[j]);
            printf("ALL PRINTED\n\n");
            MPI_Send(A[j], offset * n, MPI_DOUBLE, worker, tags[0], MPI_COMM_WORLD);
            printf("WOAH!\n\n", worker);
            MPI_Send(&B[0][i], offset * n, MPI_DOUBLE, worker, tags[1], MPI_COMM_WORLD);
            printf("send finished for worker %d \n", worker);
        }
    return 0;
}

/*
 * Description: receives results from workers
 */
double** master_receiver(int n, int offset) {
    int i, j, worker = 0;
    double** res;

    res = matrix_creator(n, n);
    for (j = 0; j < n; j += offset)
        for (i = 0; i < n; i += offset) {
            worker++;
            MPI_Recv(&res[j][i], offset * offset, MPI_DOUBLE, worker, tags[2], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }

    return res;
}

/*
 * Argv used for matrix dimension n
 */
int main(int argc, char *argv[]) {
    /*MPI variables*/
    int myrank, numnodes;

    /*Test variables
     *int ind_split, req_tag = 0, ans_tag = 2;
     *char message[100];
     */

    /*matrix variables*/
    int n = atoi(argv[1]); /*matrix n given by the user*/
    int mb, offset; /*offset is number of rows/columns for each process*/
    double** A;
    double** B;


    /*MPI initialization*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    /*variables init*/
    mb = sqrt(numnodes-1);
    offset = n / mb;

    /*show who I am*/
    printf("I'm process %d\n", myrank);

    if (myrank == 0) {
        /* matrix creation */
        A = matrix_creator(n, n);
        B = matrix_creator(n, n);

        /*init matrices with random values*/
        matrix_init(A, n);
        matrix_init(B, n);

        /*debug*/
        printf("Matices correctly created. I will print them:\n");
        printmatrix(n,n, A);
        printf("\n");
        printmatrix(n,n, B);
        printf("\n");

        /*test mpi with send message
        for (ind_split = 1; ind_split <= numnodes - 1; ind_split++) {
            sprintf(message, "This message for processor number %d\n", ind_split);
            MPI_Send(message, strlen(message) + 1, MPI_CHAR, ind_split, req_tag, MPI_COMM_WORLD);
        }*/

        /*split matrix in pieces & send matrix pieces*/
        master_sender(A, B, offset, n);

        /*debug*/
        printf("Master has just passed the master_sender funct\n");
    } else {
        /*data structure for incoming rows & cols*/
        double** rows = matrix_creator(n, offset);
        double** cols = matrix_creator(offset, n);

        /*result matrix*/
        double** res;

        /*test mpi_recv with message*/
        /*MPI_Recv(message, 100, MPI_CHAR, 0, req_tag, MPI_COMM_WORLD, &status);*/

        /*recv for rows of A and cols of B*/
        MPI_Recv(rows[0], n * offset, MPI_DOUBLE, 0, tags[0], MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(cols[0], offset * n, MPI_DOUBLE, 0, tags[1], MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /*work and free rows and cols*/
        res = mult(rows, cols, n, offset);
        freematrix(rows);
        freematrix(cols);

        /*send work back to master and free matrix*/
        MPI_Send(&res[0][0], offset * offset, MPI_DOUBLE, 0, tags[2], MPI_COMM_WORLD);
        freematrix(res);

        /*test mpi with send message
        sprintf(message, "%s COMPUTED BY PROCESSOR NUMBER: %d\n", message, myrank);
        MPI_Send(message, strlen(message) + 1, MPI_CHAR, 0, ans_tag, MPI_COMM_WORLD);
         */
    }

    MPI_Barrier(MPI_COMM_WORLD);
    /*collect all the stuff*/
    if (myrank == 0) {

        /*test mpi with send message
         * for (ind_split = 1; ind_split <= numnodes - 1; ind_split++) {
            MPI_Recv(message, 100, MPI_CHAR, ind_split, ans_tag, MPI_COMM_WORLD, &status);
            printf("Printing message: %s\n", message);
        }*/

        /*receive pieces and compute final matrix*/
        double** res = master_receiver(n, offset);

        /*print final matrix and free memory of matrices A, B and res*/
        printmatrix(n, n, res);
        freematrix(A);
        freematrix(B);
        freematrix(res);
    }

    MPI_Finalize();
    return 0;
}



