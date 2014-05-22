
/*
 * File:   mm.c
 * Author:
 *
 * Created on 12 marzo 2014, 15.20
 */

/* include header files */
#include "header.h"
#include "MpiStopwatch.h"
#include "inout.h"

/**
 * Computes multiplication
 *
 * @param rows  Pointer to the rows offsetxn matrix
 * @param cols  Pointer to the cols offsetxn matrix (it's transposed)
 * @param n     Dimension of the original matrices
 * @param offset        Dimension of the local worker's matrix
 * @return res  Pointer to the result offsetxoffset matrix
 */
double** mult(double** rows, double** cols, int n, int offset, int load) {
    int i, j, k;
    double** res; /*data structure for resulting matrix*/

    /*stopwatch start*/
    //Stopwatch watch = StopwatchCreate();
    //StopwatchStart(watch);

    res = matrix_creator(offset, offset);
    for (i = 0; i < offset; i++) {
        for (j = 0; j < offset; j++) {
            res[i][j] = 0;
            for (k = 0; k < n; k++) {
                /*notice that cols came from B that is already transposed*/
                res[i][j] += rows[i][k] * cols[j][k];
                if (load == 0)
                    res[i][j] += rows[i][k] * cols[j][k];
                else
                    res[i][j] += heavy(rows[i][k], load) * cols[j][k];
            }
        }
    }

    /*stopwatch stop*/
    /*StopwatchStop(watch);
    StopwatchPrintWithComment("Time in mult function is: %f\n", watch);
    free(watch);*/

    return res;
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
    int i, j, worker = 0;
    double * tempA, * tempB;

    /*stopwatch start*//*
    Stopwatch watch = StopwatchCreate();
    StopwatchStart(watch);*/

    for (j = 0; j < n; j += offset)
        for (i = 0; i < n; i += offset) {
            worker++; /*increment worker number*/

            /*vectorize the two pieces of matrices in order to send them*/
            tempA = matrix_vectorizer(offset, n, &A[j]);
            tempB = matrix_vectorizer(offset, n, &B[i]);

            /*MPI send of rows and cols. Notice tag 0 for rows, tag 1 for cols*/
            MPI_Send(tempA, offset * n, MPI_DOUBLE, worker, 0, MPI_COMM_WORLD);
            MPI_Send(tempB, offset * n, MPI_DOUBLE, worker, 1, MPI_COMM_WORLD);

            free(tempA);
            free(tempB);
        }

    /*stopwatch stop*/
    /*StopwatchStop(watch);
    StopwatchPrintWithComment("Time in master_sender function is: %f\n", watch);
    free(watch);*/
}

/**
 * Receives results from workers and create the final matrix
 *
 * @param n     Dimension of the final nxn matrix
 * @param offset        Dimension of one of the incoming matrices
 * @return res  Pointer to the final result nxn matrix
 */
double** master_receiver(int n, int offset) {
    int i, j, x, y, worker = 0;
    double** res, ** res_temp;

    /*stopwatch start*/
    /*Stopwatch watch = StopwatchCreate();
    StopwatchStart(watch);*/

    res = matrix_creator(n, n);
    for (j = 0; j < n; j += offset)
        for (i = 0; i < n; i += offset) {
            worker++;
            double * res_vect = (double *) malloc(offset * offset * sizeof (double));

            /*receives a piece of the final matrix in vector form, tag 2 used*/
            MPI_Recv(res_vect, offset * offset, MPI_DOUBLE, worker, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            /*transform the received vector in a matrix*/
            res_temp = devectorizer(offset, offset, res_vect);
            free(res_vect);

            /*computes the corresponding piece of the final matrix*/
            for (x = 0; x < offset; x++)
                for (y = 0; y < offset; y++)
                    res[i + x][y + j] = res_temp[x][y];
            freematrix(offset, res_temp);
        }

    /*stopwatch stop*/
    /*StopwatchStop(watch);
    StopwatchPrintWithComment("Time in master_receiver function is: %f\n", watch);
    free(watch);*/

    return res;
}

/**
 *
 * @param argc
 * @param argv Used in order to give the dimension of the matrices
 * @return 0
 */
int main(int argc, char *argv[]) {
    /*MPI variables*/
    int myrank, numnodes;

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

    // heavy load f(A) abilitation
    int load_bool = atoi(argv[3]);

    /*I\O*/
    int inout_bool = atoi(argv[2]);
    char * homePath = getenv("HOME"); /*homepath*/

    /*testfile strings*/
    char func_field[20], *io_field;
    if (load_bool == 0)
        snprintf(func_field, sizeof func_field, "low");
    else
        snprintf(func_field, sizeof func_field, "func%d", load_bool);
    if (inout_bool == 0)
        io_field = "no_io";
    else
        io_field = "io";

    /*CSV file support*/
    char *op = "mmfast", final[256];
    int myid;

    /*show who I am*/
    /*printf("I'm process %d\n", myrank);*/

    /*stopwatch start*/
    Stopwatch watch = StopwatchCreate();
    Stopwatch watchs = StopwatchCreate();

    if (myrank == 0) {

        /*stopwatch start*/
        StopwatchStart(watch);
        StopwatchStart(watchs);

        /*input type evaluation*/
        if (inout_bool == 0) { /*no input, so randomly generated matrix*/

            /* matrix creation */
            A = matrix_creator(n, n);
            B = matrix_creator(n, n);

            /*init matrices with random values*/
            simple_matrix_init(A, n);
            simple_matrix_init(B, n);

        } else if (inout_bool == 1) {
            /*input filename generation*/
            char infile[256];
            snprintf(infile, sizeof infile, "%s/hpc_temp/hpc_input/mat%d.csv", homePath, n);

            /* matrix loading */
            A = matrix_loader(infile);
            B = matrix_loader(infile);
        } else {
            printf("Error on input output value!!!\n");
            return 0;
        }


        /*transpose B for simplicity*/
        matrix_transposer(n, B);

        /*split matrix in pieces & send matrix pieces*/
        master_sender(A, B, offset, n);

        /*printf("Master has just passed the master_sender funct\n\n");*/

        /*stopwatch stop*/
        StopwatchStop(watchs);

    } else {

        /*start timer*/
        StopwatchStart(watchs);

        /*data structure for incoming rows & cols*/
        double** rows;
        double** cols;
        double * temp_rows = (double *) malloc(offset * n * sizeof (double));
        double * temp_cols = (double *) malloc(offset * n * sizeof (double));

        /*result matrix*/
        double** res;
        double * res_vect;

        /*recv for rows of A and cols of B*/
        MPI_Recv(temp_rows, offset * n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(temp_cols, n * offset, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        /*de-vectorize temp_rows & temp_cols*/
        rows = devectorizer(offset, n, temp_rows);
        cols = devectorizer(offset, n, temp_cols);

        /*work and free rows and cols*/
        res = mult(rows, cols, n, offset, load_bool);
        freematrix(offset, rows);
        freematrix(offset, cols);
        free(temp_rows);
        free(temp_cols);

        /*vectorize res in order to send it*/
        res_vect = matrix_vectorizer(offset, offset, res);

        /*send work back to master and free matrix*/
        MPI_Send(res_vect, offset * offset, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);

        /*stopwatch stop*/
        StopwatchStop(watchs);

        freematrix(offset, res);
        free(res_vect);

        /*stopwatch stop*/
        //StopwatchStop(watch);
        //StopwatchPrintWithComment("Worker total time: %f\n", watch);
    }

    //MPI_Barrier(MPI_COMM_WORLD);

    /*collect all the stuff*/
    if (myrank == 0) {
        /*start timer*/
        StopwatchStart(watchs);

        /*receive pieces and compute final matrix*/
        double** res;
        res = master_receiver(n, offset);

        /*output type evaluation*/
        if (inout_bool == 0) { /*no output, so no result (or print)*/
            /*printmatrix(n, n, res);*/
        } else {
            /*char outfile[256];
            snprintf(outfile, sizeof outfile, "%s/hpc_temp/hpc_output/%s_dim%d_nproc%d.csv", homePath, op, n, numnodes - 1);
            matrix_writer(n, res, outfile);*/
        }
        /*free memory of matrices A, B and res*/
        freematrix(n, A);
        freematrix(n, B);
        freematrix(n, res);

        /*stopwatch stop*/
        StopwatchStop(watchs);

        /*stopwatch stop*/
        StopwatchStop(watch);
        StopwatchPrintWithComment("Master total time: %f\n", watch);
        myid = (int) MPI_Wtime(); /*my id generation*/
        snprintf(final, sizeof final, "%d,%s%s,%d,%d,%s,%s", myid, op, OPTI, numnodes - 1, n, io_field, func_field); /*final string generation*/
        StopwatchPrintToFile(final, watch);

    }

    //StopwatchPrintWithComment("Slave total time: %f\n", watchs);
    myid = (int) MPI_Wtime(); /*my id generation*/
    snprintf(final, sizeof final, "%d,%s%s,%d,%d,%s,%s,%d", myid, op, OPTI, numnodes - 1, n, io_field, func_field, myrank); /*final string generation*/
    StopwatchPrintToFile2(final, watchs);

    free(watch);
    free(watchs);
    MPI_Finalize();
    return 0;
}