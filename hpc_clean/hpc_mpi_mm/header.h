/* 
 * File:   functions.h
 * Author: jian
 *
 * Created on 13 marzo 2014, 11.42
 */

#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

/*includes*/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define TAG 13
const int tags[3] = {1, 2, 3};


/*
 * 
 * FUNCTION DEFINITION
 */

/* 
 * Creates axb matrix
 */
double** matrix_creator(int a, int b) {
    double **mat = (double **) malloc(a * sizeof (double*));
    int i;
    for (i = 0; i < a; i++)
        mat[i] = (double *) malloc(b * sizeof (double));
    return mat;
}

int matrix_init(double** mat, int n) {
    int i, j; /*matrix indexes*/

    /* random seed initialization */
    srand(time(NULL));

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            mat[i][j] = ((double) rand() / (double) RAND_MAX);
    }
    return 0;
}

/*
 Initialize a matrix with simple stram of integers (casted in double)
 */
int simple_matrix_init(double** mat, int n) {
    int i, j; /*matrix indexes*/
    int counter=1;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            mat[i][j] = counter++;
    }
    return 0;
}

/*
 * free a matrix
 */
void freematrix(int n, double** mat) {
    int i;
    for (i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

/*print matrix axb*/
void printmatrix(int a, int b, double** C) {
    int i, j;
    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++)
            printf("%f ", C[i][j]);
        printf("\n");
    }
    printf("\n");
}

/*print vector */
void printvector(int a, double* C) {
    int i, j;
    for (i = 0; i < a; i++) 
        printf("%f ", C[i]);
    printf("\n");
}

/*
 transform a matrix in a vector
 */
double * matrix_vectorizer(int a, int b, double ** mat){
    double * vet = (double *) malloc(a * b * sizeof (double));
    int i, j, offset = 0;
    for(i=0;i<a;i++){
        for(j=0;j<b;j++){
            vet[j+offset] = mat[i][j];
        }
        offset = offset + b;
    }
    return vet;
}

/*
 transform a vector in a matrix
 */
double ** devectorizer(int a, int b, double * vet){
    double ** mat = matrix_creator(a,b);
    int i, j, offset = 0;
    
    for(i=0;i<a;i++){
        for(j=0;j<b;j++){
            mat[i][j] =  vet[j+offset];
        }
    offset = offset + b;
    }
    
    return mat;
}

/*
 transpose a square matrix
 */
int matrix_transposer(int n, double ** A) {
    int i, j; /*indexes*/
    double temp;
    for (i = 0; i < n-1; i++)
        for (j = i + 1; j < n; j++) {
            temp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = temp;
        }
    return 0;
}

#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}
#endif

#endif	/* FUNCTIONS_H */

