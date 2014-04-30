/* 
 * File:   matcreator.c
 * Author: 
 *
 * Created on 29 aprile 2014, 17.35
 */

/* include header files */
#include "header.h"
#include "inout.h"

/*
 * 
 */
int main(int argc, char** argv) {
    double ** mat;
    int dim;
    char filename[256];

    /*getting dimension from argv[1]*/
    dim = atoi(argv[1]);

    /*matrix creation*/
    mat = matrix_creator(dim, dim);
    if (mat == NULL) {
        printf("Error in matrix creation!!\n");
        return 0;
    }
    
    /*matrix initialization*/
    matrix_init(mat, dim);
    
    /*filename creation*/
    snprintf(filename, sizeof filename, "mat%d.csv", dim);

    /*matrix saving*/
    if (matrix_writer(dim, mat, filename) == 0) {
        printf("Error in matrix writing!!\n");
        return 0;
    }
    
    printf("A file .csv with dimension %dx%d with name %s has successfully been created!!\n", dim, dim, filename);
    return (EXIT_SUCCESS);
}