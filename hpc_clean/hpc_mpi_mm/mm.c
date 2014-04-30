/* 
 * File:   cannon.c
 * Author: asus
 *
 * Created on 14 aprile 2014, 10.56
 */

#include "header.h"
#include "MpiStopwatch.h"
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

    //freematrix(n, r_swap);
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

    //freematrix(n, c_swap);
}

double ** matrix_block(double** matrix, int n, int nblock) {
    int i, j, k, x, y, offset = n / sqrt(nblock);
    double * tempM, **block, ** Mblock;

    block = matrix_creator(offset, offset);
    Mblock = matrix_creator(nblock, (offset * offset));

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

            //free(tempM);
        }

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

void matrix_mult(double** A, double** B, double** C, int dim) {
    int i, j, k, l;

    for (i = 0; i < dim; i++) {
        l = 0;
        for (j = 0; j < dim; j++) {
            for (k = 0; k < dim; k++) {
                C[i][j] += A[i][k] * B[l][k];
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

    double **A, **B, **C, *tmpA, *tmpB, **Ablock, **Bblock;
    double startTime, endTime;
    int master, nblock, numElements, offset, myrank, rankR, rankC, numnodes, N, i, j, k, l;
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

    // allocate A, B, and C --- note that you want these to be
    // contiguously allocated.  Workers need less memory allocated.

    if (myrank == master) {
        printf("Printf atoi N: %d\n", N);
        printf("Printf numnodes: %d\n", numnodes);

        /*start timer*/
        StopwatchStart(watch);

        A = matrix_creator(N, N);

        B = matrix_creator(N, N);

        C = matrix_creator(N, N);

        printf("MASTER. A,B,C allocated\n");

        // initialize A and B
        simple_matrix_init(A, N);
        simple_matrix_init(B, N);
        zero_matrix_init(C, N, N);

        // suddivisione in blocchi della matrice
        Ablock = matrix_block(A, N, nblock);
        Bblock = matrix_block(B, N, nblock);

        //printmatrix(nblock, numElements, Ablock);
        //printmatrix(nblock, numElements, Bblock);

        numElements = dim*dim;

        // send
        for (i = 0; i < nblock; i++) {
            MPI_Send(Ablock[i], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
            MPI_Send(Bblock[i], numElements, MPI_DOUBLE, i, TAG, MPI_COMM_WORLD);
        }

        printf("MASTER. Pieces of A and B sent.\n");

    } else {
        // receive my part of A and B
        numElements = dim*dim;
        lato = (int) sqrt(nblock);

        Ablock = matrix_creator(lato, numElements);
        Bblock = matrix_creator(lato, numElements);

        double **Aswap = matrix_creator(1, numElements);
        double **Bswap = matrix_creator(1, numElements);

        MPI_Recv(Aswap[0], numElements, MPI_DOUBLE, master, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(Bswap[0], numElements, MPI_DOUBLE, master, TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        printf("Myrank is %d. Pieces of A and B received.\n", myrank);

        C = matrix_creator(dim, dim);

        // initialize C
        zero_matrix_init(C, dim, dim);

        for (i = 0; i < numElements; i++) {
            Ablock[0][i] = Aswap[0][i];
            Bblock[0][i] = Bswap[0][i];
        }

        for (i = 1; i < lato; i++) {

            row_dest = getRankRowDest(myrank, nblock);
            row_mit = getRankRowMit(myrank, nblock);

            // invio e ricezione del blocco A
            MPI_Sendrecv_replace(Aswap[0], numElements, MPI_DOUBLE, row_dest, 1, row_mit, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            col_dest = getRankColDest(myrank, nblock);
            col_mit = getRankColMit(myrank, nblock);

            // invio e ricezione del blocco B
            MPI_Sendrecv_replace(Bswap[0], numElements, MPI_DOUBLE, col_dest, 1, col_mit, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            for (j = 0; j < numElements; j++) {
                Ablock[i][j] = Aswap[0][j];
                Bblock[i][j] = Bswap[0][j];
            }
        }


        /*

        for (index = 0; index < lato; index++) {
            A = devectorizer(dim, dim, Ablock[0]);
            B = devectorizer(dim, dim, Bblock[0]);
            matrix_transposer(dim, B);

            // Multiplication
            matrix_mult(A, B, C, dim);

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

        //stopwatch stop
        StopwatchStop(watch);
        StopwatchPrintWithComment("Master total time: %f\n\n", watch);

        printf("MASTER. Print C.\n");
        printmatrix(N, N, C);*/
    }

    MPI_Finalize();

    return 0;
}
