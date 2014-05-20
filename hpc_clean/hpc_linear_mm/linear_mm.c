/*
 * File:   cannon.c
 * Author: asus
 *
 * Created on 14 aprile 2014, 10.56
 */

#include "header.h"
#include "MpiStopwatch.h"
#include "inout.h"
#define TAG 13

void zero_matrix_init(double** mat, int a, int b) {
    int i, j; /*matrix indexes*/

    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++)
            mat[i][j] = 0.0;
    }
}

void matrix_mult(double** A, double** B, double** C, int r, int c, int load) {
    int i, j, k, l;

    for (i = 0; i < r; i++) {
        l = 0;
        for (j = 0; j < r; j++) {
            for (k = 0; k < c; k++) {
                if (load == 0)
                    C[i][j] += A[i][k] * B[l][k];
                else
                    C[i][j] += heavy(A[i][k], load) * B[l][k];
            }
            l++;
        }
    }
}

int main(int argc, char** argv) {

    double **A, **B, **C;
    int numnodes, N;

    /*stopwatch*/
    Stopwatch watch = StopwatchCreate();

    N = atoi(argv[1]);
    MPI_Init(&argc, &argv);

    /*I\O*/
    int inout_bool = atoi(argv[2]);
    char * homePath = getenv("HOME"); /*homepath*/

    // heavy load f(A) abilitation
    int load_bool = atoi(argv[3]);

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
    char *op = "l_mm", final[256];
    int myid;

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

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

    matrix_transposer(N, B);

    C = matrix_creator(N, N);

    printf("MASTER. A,B,C allocated\n");

    zero_matrix_init(C, N, N);


    // Multiplication
    matrix_mult(A, B, C, N, N, load_bool);


    /*output type evaluation*/
    if (inout_bool == 0) { /*no output, so no result (or print)*/
        /*printmatrix(N, N, res);*/
    } else {
        /*char outfile[256];
        snprintf(outfile, sizeof outfile, "%s/hpc_temp/hpc_output/%s_dim%d_nproc%d.csv", homePath, op, N, numnodes - 1);
        matrix_writer(N, C, outfile);*/
    }

    /*stopwatch stop*/
    StopwatchStop(watch);
    StopwatchPrintWithComment("Master total time: %f\n", watch);
    myid = (int) MPI_Wtime(); /*my id generation*/
    snprintf(final, sizeof final, "%d,%s%s,%d,%d,%s,%s", myid, op, OPTI, numnodes - 1, N, io_field, func_field); /*final string generation*/
    StopwatchPrintToFile(final, watch);

    // free
    freematrix(N, A);
    freematrix(N, B);
    freematrix(N, C);
    free(watch);
    
    MPI_Finalize();

    return 0;
}
