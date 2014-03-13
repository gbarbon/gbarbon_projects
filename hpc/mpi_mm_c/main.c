/* 
 * File:   main.c
 * Author: jian
 *
 * Created on 13 marzo 2014, 15.43
 */

#include "functions.h"   

/*
 * 
 */
int main(int argc, char** argv) {
    int N = atoi(argv[1]);
    
    MPI_Init(&argc, &argv);
    mm(N);
    
    
    return (EXIT_SUCCESS);
}

