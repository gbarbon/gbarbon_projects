#include "header.h"

int* coordinate(int procNum, int totalProc) {
    int* coord = (int*) calloc(2, sizeof (int)); //aggiunto (int*)
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    return coord;
}

void printmatrix(int N, double** C) {
    // print matrix

    int i, j; //matrix indexes

    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            printf("%f ", C[i][j]);
        printf("\n");

    }
}
