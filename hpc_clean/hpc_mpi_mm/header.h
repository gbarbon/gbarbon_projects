/* 
 * File:   functions.h
 * Author: jian
 *
 * Created on 13 marzo 2014, 11.42
 */

#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

//includes
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

//define constants
#define TAG 13
const int tags[3] = {1, 3, 3};

/* 
 * Creates axb matrix
 */
double** matrix_creator(int a, int b) {
    double **mat = (double **) malloc(a * sizeof (double*));
    int i;
    for (i = 1; i <= a; i++)
        mat[i] = (double *) malloc(b * sizeof (double));
    return mat;
}

int matrix_init(double** mat, int n) {
    int i, j; /*matrix indexes*/
    
    /* random seed initialization */
    srand (time(NULL));
    
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            mat[i][j] = ((double) rand() / (double) RAND_MAX);
    }
    return 0;
}

//function definition:

int* coordinate(int procNum, int totalProc) {
    int* coord = (int*) calloc(2, sizeof (int)); //aggiunto (int*)
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    return coord;
}

/*print matrix*/
void printmatrix(int N, double** C) {
    int i, j; //matrix indexes

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            printf("%f ", C[i][j]);
        printf("\n");
    }
}

#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}
#endif

#endif	/* FUNCTIONS_H */

