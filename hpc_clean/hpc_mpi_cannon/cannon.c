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

void skewing_row(double ** M, int n) {
    int i, j, k, index;

    double **r_swap = (double **) malloc(sizeof (double*) * n);
    for (i = 0; i < n; i++) {
        r_swap[i] = M[i];
    }

    k = 1;
    for (i = sqrt(n); i < n; i = i + sqrt(n)) {
        for (j = 0; j < sqrt(n); j++) {
            index = (j + k) % (int) sqrt(n);
            M[i + j] = r_swap[i + index];
        }
        k++;
    }
}

void skewing_column(double ** M, int n) {
    int i, j, k, index;

    double **c_swap = (double **) malloc(sizeof (double*) * n);
    for (i = 0; i < n; i++) {
        c_swap[i] = M[i];
    }

    k = 1;
    for (i = 1; i < sqrt(n); i++) {
        if (i % (int) sqrt(n) != 0) {
            for (j = 0; j < sqrt(n); j++) {
                index = (j + k) % (int) sqrt(n);
                M[i + (int) sqrt(n) * j] = c_swap[i + (int) sqrt(n) * index];
            }
            k++;
        }
    }
}

double ** matrix_block(double** matrix, int n, int nblock) {
    int i, j, k, x, y, offset = n / sqrt(nblock);
    double * tempM, **block, ** Mblock;

    block = matrix_creator(offset, offset);
    Mblock = (double **) malloc(nblock * sizeof (double*));

    k = 0;
    for (i = 0; i < n; i += offset)
        for (j = 0; j < n; j += offset) {

            for (x = i; x < offset + i; x++)
                for (y = j; y < offset + j; y++) {
                    block[x - i][y - j] = matrix[x][y];
                }

            /*vectorize the two pieces of matrices*/
            tempM = matrix_vectorizer(offset, offset, block);
            Mblock[k] = tempM;
            k++;
        }

    freematrix(offset, block);

    return Mblock;
}

void block_matrix(double ** matrix, double* vett, int block, int n) {
    int i, j, x, y, el = 0, dim = n / sqrt(block);

    /*base point (in final matrix) scrolling*/
    /*i.e. 0,0 - 0,2 - 2,0 - 2,2 with N=4 and nproc=4 */
    for (i = 0; i < n; i += dim)
        for (j = 0; j < n; j += dim)
            /*block scrolling, dim: block dimension*/
            for (x = i; x < dim + i; x++)
                for (y = j; y < dim + j; y++) {
                    //printf("x: %d, y: %d, el: %f\n", x, y, C_vett[el]);
                    matrix[x][y] = vett[el];
                    el++;
                }
}

void zero_matrix_init(double** mat, int a, int b) {
    int i, j; /*matrix indexes*/

    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++)
            mat[i][j] = 0.0;
    }
}

void matrix_mult(double** A, double** B, double** C, int dim, int load) {
    int i, j, k, l;

    for (i = 0; i < dim; i++) {
        l = 0;
        for (j = 0; j < dim; j++) {
            for (k = 0; k < dim; k++) {
                if (load == 0)
                    C[i][j] += A[i][k] * B[l][k];
                else
                    C[i][j] += heavy(A[i][k], load) * B[l][k];
            }
            l++;
        }
    }
}

// determinazione del rank del processo a cui devo inviare il sottoblocco di A posseduto

int getRankRowDest(int rank, int np) {
    int rank_dest;
    int proc_lato = (int) sqrt(np);
    // se siamo ad inizio riga
    if ((rank % proc_lato) == 0)
        rank_dest = rank + (proc_lato - 1);
    else
        rank_dest = rank - 1;

    return rank_dest;
}

// determinazione del rank del processo da cui devo ricevere il sottoblocco di A

int getRankRowMit(int rank, int np) {
    int rank_mit;
    int proc_lato = (int) sqrt(np);
    // se siamo sull'ultima colonna
    if ((rank + 1) % proc_lato == 0)
        rank_mit = rank - (proc_lato - 1);
    else
        rank_mit = rank + 1;

    return rank_mit;
}

// determinazione del rank del processo a cui devo inviare il sottoblocco di B

int getRankColDest(int rank, int np) {
    int rank_dest;
    int proc_lato = (int) sqrt(np);

    if (rank < proc_lato)
        rank_dest = np - (proc_lato - rank);
    else
        rank_dest = rank - proc_lato;

    return rank_dest;
}

// determinazione del rank del processo da cui ricevere il sottoblocco di B

int getRankColMit(int rank, int np) {
    int rank_mit;
    int proc_lato = (int) sqrt(np);

    if (rank >= (np - proc_lato))
        rank_mit = (rank + proc_lato) % np;
    else
        rank_mit = rank + proc_lato;

    return rank_mit;
}

int main(int argc, char** argv) {

    double **A, **B, **C, **Ablock, **Bblock;
    int master, nblock, numElements, offset, myrank, numnodes, N, i;
    int row_dest, row_mit, col_dest, col_mit, index, lato, dim;

    /*stopwatch*/
    Stopwatch watch = StopwatchCreate();

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &numnodes);

    N = atoi(argv[1]);
    nblock = numnodes - 1;
    master = nblock;

    dim = N / sqrt(nblock);

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
    char *op = "cannon", final[256];
    int myid;

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == master) {
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

        printf("MASTER. A,B,C allocated\n");

        zero_matrix_init(C, N, N);

        // suddivisione in blocchi della matrice
        Ablock = matrix_block(A, N, nblock);
        Bblock = matrix_block(B, N, nblock);

        numElements = dim*dim;

        printmatrix(nblock, numElements, Ablock);
        printmatrix(nblock, numElements, Bblock);

        // skewing        
        skewing_row(Ablock, nblock);
        skewing_column(Bblock, nblock);

        printmatrix(nblock, numElements, Ablock);
        printmatrix(nblock, numElements, Bblock);

        // send
        for (i = 0; i < nblock; i++) {
            MPI_Send(Ablock[i], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(Bblock[i], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
        }

        printf("MASTER. Pieces of A and B sent.\n");

    } else {
        // receive my part of A and B
        numElements = dim*dim;
        Ablock = matrix_creator(1, numElements);
        Bblock = matrix_creator(1, numElements);

        MPI_Recv(Ablock[0], numElements, MPI_DOUBLE, master, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(Bblock[0], numElements, MPI_DOUBLE, master, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Myrank is %d. Pieces of A and B received.\n", myrank);

        C = matrix_creator(dim, dim);

        // initialize C
        zero_matrix_init(C, dim, dim);

        lato = (int) sqrt(nblock);

        for (index = 0; index < lato; index++) {
            A = devectorizer(dim, dim, Ablock[0]);
            B = devectorizer(dim, dim, Bblock[0]);
            matrix_transposer(dim, B);

            // Multiplication
            matrix_mult(A, B, C, dim, load_bool);

            row_dest = getRankRowDest(myrank, nblock);
            row_mit = getRankRowMit(myrank, nblock);

            // invio e ricezione del blocco A
            MPI_Sendrecv_replace(Ablock[0], numElements, MPI_DOUBLE, row_dest, 1, row_mit, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            col_dest = getRankColDest(myrank, nblock);
            col_mit = getRankColMit(myrank, nblock);

            // invio e ricezione del blocco B
            MPI_Sendrecv_replace(Bblock[0], numElements, MPI_DOUBLE, col_dest, 1, col_mit, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        }

        double *C_vett = matrix_vectorizer(dim, dim, C);
        MPI_Send(C_vett, numElements, MPI_DOUBLE, master, TAG, MPI_COMM_WORLD);

        // free
        freematrix(dim, A);
        freematrix(dim, B);
        freematrix(dim, C);
        freematrix(1, Ablock);
        freematrix(1, Bblock);
        free(C_vett);
    }

    if (myrank == master) {
        offset = 0;

        double **tempC = matrix_creator(nblock, numElements);

        for (i = 0; i < nblock; i++) {
            MPI_Recv(tempC[offset], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            offset++;
        }

        double *C_vett = matrix_vectorizer(nblock, numElements, tempC);

        block_matrix(C, C_vett, nblock, N);

        /*output type evaluation*/
        if (inout_bool == 0) { /*no output, so no result (or print)*/
            /*printmatrix(n, n, res);*/
        } else {
            /*char outfile[256];
            snprintf(outfile, sizeof outfile, "%s/hpc_temp/hpc_output/%s_dim%d_nproc%d.csv", homePath, op, N, numnodes - 1);
            matrix_writer(N, C, outfile);*/
        }

        /*stopwatch stop*/
        StopwatchStop(watch);
        StopwatchPrintWithComment("Master total time: %f\n", watch);
        myid = (int) MPI_Wtime(); /*my id generation*/
        snprintf(final, sizeof final, "%d,%s%s,%d,%d,%s,%s", myid, op, OPTI, numnodes - 1, N,  io_field, func_field); /*final string generation*/
        StopwatchPrintToFile(final, watch);

        // free
        freematrix(N, A);
        freematrix(N, B);
        freematrix(N, C);
        freematrix(nblock, Ablock);
        freematrix(nblock, Bblock);
        freematrix(nblock, tempC);
        free(C_vett);
    }

    free(watch);
    MPI_Finalize();

    return 0;
}
