

void  printmatrix() {
    // print matrix
    
    int i, j; //matrix indexes
    
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++)
            printf("%f ", C[i][j]);
        printf("\n");

    }

