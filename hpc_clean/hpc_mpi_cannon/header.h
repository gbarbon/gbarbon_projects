/* 
 * File:   functions.h
 * Author: 
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

/* CONSTANTS DEFINITION */
/*const int tags[3] = {1, 2, 3};*/

/* FUNCTION DEFINITION */

/**
 * Creates axb matrix
 * 
 * @param a     Number of rows of the matrix
 * @param b     Number of columns of the matrix
 * @return mat  Pointer to a double axb matrix
 */
double** matrix_creator(int a, int b) {
    double **mat = (double **) malloc(a * sizeof (double*));
    int i;
    for (i = 0; i < a; i++)
        mat[i] = (double *) malloc(b * sizeof (double));
    return mat;
}

/**
 * Initializes a matrix with random values
 * 
 * @param mat   Pointer to a nxn square matrix
 * @param n     Dimension of the matrix
 */
void matrix_init(double** mat, int n) {
    int i, j; /*matrix indexes*/

    srand(time(NULL)); /* random seed initialization */

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            mat[i][j] = ((double) rand() / (double) RAND_MAX);
    }
}

/**
 * Initializes a matrix with simple stream of integers (casted in double)
 * 
 * @param mat   Pointer to a nxn square matrix
 * @param n     Dimension of the matrix
 */
void simple_matrix_init(double** mat, int n) {
    int i, j; /*matrix indexes*/
    int counter = 1;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            mat[i][j] = counter++;
    }
}

/**
 * Frees a matrix
 * 
 * @param n     Rows of the matrix
 * @param mat   Pointer to a matrix
 */
void freematrix(int n, double** mat) {
    int i;
    for (i = 0; i < n; i++) {
        free(mat[i]);
    }
    free(mat);
}

/**
 * Prints matrix axb
 * 
 * @param a     Number of rows of the matrix
 * @param b     Number of cols of the matrix
 * @param C     Pointer to a axb matrix
 */
void printmatrix(int a, int b, double** C) {
    if (a <= 10 && b <= 10) {
        int i, j;
        for (i = 0; i < a; i++) {
            for (j = 0; j < b; j++)
                printf("%f ", C[i][j]);
            printf("\n");
        }
        printf("\n");
    }
}

/**
 * Prints a vector
 * 
 * @param a     Vector dimension
 * @param C     Pointer to a vector
 */
void printvector(int a, double* C) {
    int i, j;
    for (i = 0; i < a; i++)
        printf("%f ", C[i]);
    printf("\n");
}

/**
 * Transforms a matrix in a vector
 * 
 * @param a     Number of rows of the matrix    
 * @param b     Number of cols of the matrix
 * @param mat   Pointer to the axb matrix
 * @return vet  Vector composed of all the rows of the matrix
 */
double * matrix_vectorizer(int a, int b, double ** mat) {
    double * vet = (double *) malloc(a * b * sizeof (double));
    int i, j, offset = 0;
    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            vet[j + offset] = mat[i][j];
        }
        offset = offset + b;
    }
    return vet;
}

/**
 * Transforms a vector in a matrix
 * 
 * @param a     Number of rows of the matrix
 * @param b     Number of cols of the matrix
 * @param vet   Pointer to the vector
 * @return mat  Matrix axb
 */
double ** devectorizer(int a, int b, double * vet) {
    double ** mat = matrix_creator(a, b);
    int i, j, offset = 0;

    for (i = 0; i < a; i++) {
        for (j = 0; j < b; j++) {
            mat[i][j] = vet[j + offset];
        }
        offset = offset + b;
    }

    return mat;
}

/**
 * Transpose a square matrix
 * 
 * @param n     Dimension of the matrix
 * @param A     Pointer to the nxn matrix
 */
void matrix_transposer(int n, double ** A) {
    int i, j; /*indexes*/
    double temp;
    for (i = 0; i < n - 1; i++)
        for (j = i + 1; j < n; j++) {
            temp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = temp;
        }
}

double heavy(double a) {
    int i;
    for (i = 0; i < 5; i++)
        a = pow(a, i);
    return a;
}

#endif	/* FUNCTIONS_H */
