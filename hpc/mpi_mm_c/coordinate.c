int* coordinate(int procNum, int totalProc) {
    int* coord = (int*) calloc(2, sizeof (int)); //aggiunto (int*)
    int var;
    var = sqrt(totalProc);
    coord[0] = procNum / var;
    coord[1] = procNum % var;
    return coord;
}
