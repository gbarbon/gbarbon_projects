/* 
 * File:  
 * Author: 
 *
 * 
 */

#ifndef INOUT_H
#define	INOUT_H

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
            mat[i][j++] = atof(record);
            record = strtok(NULL, ",");
        }
        i++;
    }
    
    return mat;
}

#endif	/* INOUT_H */