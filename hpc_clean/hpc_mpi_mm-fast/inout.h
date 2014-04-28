/* 
 * File:  
 * Author: 
 *
 * 
 */

#ifndef INOUT_H
#define	INOUT_H

/*includes*/
#include "header.h"

/* FUNCTION DEFINITION */

/**
 * Loads axb matrix
 * 
 * @param input_file     Name of the input file
 * @return mat  Pointer to a double nxn matrix
 */
double ** matrix_loader(char* input_file) {
    char buffer[1024];
    char *record, *line;
    int i = 0, j = 0, already_col = 0;

    /*first file loading*/
    FILE *inputf = fopen(input_file, "r");
    if (inputf == NULL) {
        printf("Impossible to load input files\n");
        return NULL;
    }

    /*this while cycle determines the dimensions of the matrix*/
    while ((line = fgets(buffer, sizeof (buffer), inputf)) != NULL) {
        if (!already_col) {
            record = strtok(line, ",");
            /*already_row used to calculate only one time the row dimension*/
            while (record != NULL) {
                j++;
                record = strtok(NULL, ",");
            }
            already_col = 1;
        }
        i++;
    }

    fclose(inputf);
    if (i != j) {
        printf("input matrix is not sqared. Return NULL\n");
        return NULL;
    }

    /*second file loading*/
    double ** mat = matrix_creator(i, j);
    i = 0;
    j = 0;
    inputf = fopen(input_file, "r");

    /*this while cycle loads the matrix*/
    while ((line = fgets(buffer, sizeof (buffer), inputf)) != NULL) {
        record = strtok(line, ",");
        while (record != NULL) {
            printf("record : %s ", record); //here you can put the record into the array as per your requirement.
            mat[i][j] = atof(record);
            //printf("Matrix %f ", mat[i][j]);
            record = strtok(NULL, ",");
            j++;
        }
        j = 0;
        i++;
    }
    fclose(inputf);
    return mat;
}

/**
 * Write dimxdim matrix into a file
 * 
 * @param dim   Matrix dimension 
 * @param mat   Pointer to a double nxn matrix
 * @param output_file   Name of the destination file
 * @return 0 on failure, 1 on success
 */
int matrix_writer(int dim, double ** mat, char* output_file) {
    int i, j;
    char buf[256];

    FILE* outf=fopen(output_file, "w");
    if (outf == NULL) {
        printf("Impossible to create output file\n");
        return 0;
    }
    
    for (i = 0; i < dim; i++)
        for (j = 0; j < dim; j++) {
            if (j < dim - 1) {
//                snprintf(buf, sizeof buf , "%f%s", mat[i][j], ",");
                snprintf(buf, sizeof buf , "%f,", mat[i][j]);
            } else {
//                snprintf(buf, sizeof buf , "%f%s", mat[i][j], "\n");
                snprintf(buf, sizeof buf , "%f\n", mat[i][j]);
            }
            fputs(buf, outf);
        }
    
    fclose(outf);
    return 1;
}


#endif	/* INOUT_H */