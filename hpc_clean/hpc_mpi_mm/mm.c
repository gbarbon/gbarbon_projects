
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
    for (i = 0; i < offset; i++) {
        for (j = 0; j < offset; j++) {
            res[i][j] = 0;
            for (k = 0; k < n; k++) {
                res[i][j] += rows[i][k] * cols[j][k];
//                printf("I: %d, J: %d, ", i, j);
//                printf("rows: %f, ", rows[i][k]);
//                printf("cols: %f, ", cols[j][k]);
//                printf("Res: %f\n", res[i][j]);
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
    double * tempA, * tempB;
    for (j = 0; j < n; j += offset)
        for (i = 0; i < n; i += offset) {
            worker++;
            //printf("\n\n");
            //            if (worker == 1) {
            //                printf("FOR node0%d: Print part of matrix A\n", worker);
            //                printmatrix(offset, n, &A[j]);
            //                printf("FOR node0%d: Print part of matrix B\n", worker);
            //                printmatrix(n, offset, &(&B[0])[i]);
            //            }
            //printf("Cycle for worker %d , offset is: %d, j now is: %d, i now is: %d \n", worker, offset, j, i);
            tempA = matrix_vectorizer(offset, n, &A[j]);
            //printf("A VECTORIZED for worker: %d\n", worker);
            tempB = matrix_vectorizer(offset, n, &B[i]);
            //            if (worker == 1) {
            //                printf("FOR node0%d: Print part of temp A\n", worker);
            //                printvector(offset*n, tempA);
            //                printf("FOR node0%d: Print part of temp B\n", worker);
            //                printvector(n*offset, tempB);
            //            }
            //printf("\n");


            //printf("B VECTORIZED for worker: %d\n", worker);

            MPI_Send(tempA, offset * n, MPI_DOUBLE, worker, 0, MPI_COMM_WORLD);
            //printf("FOR node0%d: WOAH!\n", worker);
            MPI_Send(tempB, offset * n, MPI_DOUBLE, worker, 1, MPI_COMM_WORLD);
            //printf("FOR node0%d: both send finished\n", worker);
            free(tempA);
            free(tempB);
            //printf("temp freed \n");
        }
    return 0;
}

/*
 * Description: receives results from workers
 */
double** master_receiver(int n, int offset) {
    int i, j, x, y, worker = 0;
    double** res, ** res_temp;

    res = matrix_creator(n, n);
    for (j = 0; j < n; j += offset)
        for (i = 0; i < n; i += offset) {
            worker++;
            double * res_vect = (double *) malloc(offset * offset * sizeof (double));

            MPI_Recv(res_vect, offset * offset, MPI_DOUBLE, worker, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            res_temp = devectorizer(offset, offset, res_vect);
            free(res_vect);

            for (x = 0; x < offset; x++)
                for (y = 0; y < offset; y++)
                    res[i + x][y + j] = res_temp[x][y];
            freematrix(offset, res_temp);
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
    double **A;
    double **B;


    /*MPI initialization*/
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    /*variables init*/
    mb = sqrt(numnodes - 1);
    offset = n / mb;

    /*show who I am*/
    //printf("I'm process %d\n", myrank);

    if (myrank == 0) {
        /* matrix creation */
        A = matrix_creator(n, n);
        B = matrix_creator(n, n);

        /*init matrices with random values*/
        simple_matrix_init(A, n);
        simple_matrix_init(B, n);

        /*transpose B for simplicity*/
        matrix_transposer(n, B);

        /*debug*/
        printf("Matices correctly created. I will print them:\n");
        printmatrix(n, n, A);
        printf("\n");
        printmatrix(n, n, B);
        printf("\n");

        /*test mpi with send message
        for (ind_split = 1; ind_split <= numnodes - 1; ind_split++) {
            sprintf(message, "This message for processor number %d\n", ind_split);
            MPI_Send(message, strlen(message) + 1, MPI_CHAR, ind_split, req_tag, MPI_COMM_WORLD);
        }*/

        /*split matrix in pieces & send matrix pieces*/
        master_sender(A, B, offset, n);

        /*debug*/
        printf("Master has just passed the master_sender funct\n\n");
    } else {
        /*data structure for incoming rows & cols*/
        double** rows;
        double** cols;
        double * temp_rows = (double *) malloc(offset * n * sizeof (double));
        double * temp_cols = (double *) malloc(offset * n * sizeof (double));
        /*result matrix*/
        double** res;
        double * res_vect;

        /*test mpi_recv with message*/
        //        /*MPI_Recv(message, 100, MPI_CHAR, 0, req_tag, MPI_COMM_WORLD, &status);*/
        //        if (myrank == 1) {
        //            printf("This is offset: %d\n", offset);
        //            printf("node0%d: Printing empty rows\n", myrank);
        //            printmatrix(offset, n, rows);
        //            printf("node0%d: Printing empty cols\n", myrank);
        //            printmatrix(n, offset, cols);
        //            printf("BLABLA\n\n");
        //        }
        /*recv for rows of A and cols of B*/

        //         printf("PRINT BEFORE RECV\n", myrank);
        //        if (myrank == 1) {
        //            printf("This is offset: %d\n", offset);
        //            printf("node0%d: Printing temp rows EMPTY\n", myrank);
        //            printvector(offset*n, temp_rows);
        //            printf("node0%d: Printing temp cols EMPTY\n", myrank);
        //            printvector(n*offset, temp_cols);
        //            printf("____\n\n");
        //        }

        //printf("MY RANK IS: %d\n", myrank);
        //printf("node0%d: Waiting for incoming rows\n", myrank);
        //printf("Offset is: %d , n is %d \n", offset, n);
        MPI_Recv(temp_rows, offset * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("node0%d: Waiting for incoming cols\n", myrank);
        MPI_Recv(temp_cols, n * offset, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("node0%d: Rows & cols received\n", myrank);

        //int number;
        //MPI_Recv(&number, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        //        printf("MY RANK IS: %d\n", myrank);
        //        //printf("My number is: %d\n", number);
        //        
        //            printf("This is offset: %d\n", offset);
        //            printf("node0%d: Printing temp rows\n", myrank);
        //            printvector(offset*n, temp_rows);
        //            printf("node0%d: Printing temp cols\n", myrank);
        //            printvector(n*offset, temp_cols);
        //            printf("BLABLA\n\n");

        //MPI_Barrier(MPI_COMM_WORLD);

        /*de-vectorize temp_rows & temp_cols*/
        rows = devectorizer(offset, n, temp_rows);
        cols = devectorizer(offset, n, temp_cols);

//        printf("This is offset: %d\n", offset);
//        printf("node0%d: Printing empty rows\n", myrank);
//        printmatrix(offset, n, rows);
//        printf("node0%d: Printing empty cols\n", myrank);
//        printmatrix(offset, n, cols);
//        printf("\n\n");

        /*work and free rows and cols*/
        res = mult(rows, cols, n, offset);
        freematrix(offset, rows);
        freematrix(offset, cols);
        free(temp_rows);
        free(temp_cols);


        /*vectorize res in order to send it*/
        res_vect = matrix_vectorizer(offset, offset, res);

        /*send work back to master and free matrix*/
        MPI_Send(res_vect, offset * offset, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
        freematrix(offset, res);
        free(res_vect);

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
        double** res;
        res = master_receiver(n, offset);

        /*print final matrix and free memory of matrices A, B and res*/
        printmatrix(n, n, res);
        freematrix(n, A);
        freematrix(n, B);
        freematrix(n, res);
    }

    MPI_Finalize();
    return 0;
}



