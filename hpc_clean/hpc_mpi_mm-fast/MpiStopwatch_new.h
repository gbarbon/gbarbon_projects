/*
 * Stopwatch.h
 *
 *  Created on: 
 *      Author: 
 */

#ifndef STOPWATCH_H_
#define STOPWATCH_H_

#include <mpi.h>

typedef struct {
    double _start;
    double _end;
} _stopwatch;

typedef _stopwatch* Stopwatch;

Stopwatch StopwatchCreate() {
    Stopwatch res = (Stopwatch) malloc(sizeof (_stopwatch));
    res->_start = 0;
    res->_end = 0;
    return res;
}

void StopwatchStart(Stopwatch chrono) {
    chrono->_start = MPI_Wtime(); //clock(); //time(0);
}

void StopwatchStop(Stopwatch chrono) {
    chrono->_end = MPI_Wtime();
}

float StopwatchGetElapsed(Stopwatch chrono) {
    return ((float) (chrono->_end - chrono->_start));
}

int StopwatchIsValid(Stopwatch chrono) {
    return chrono->_end >= chrono->_start;
}

void StopwatchPrintResolution() {

    double tick = MPI_Wtick() * 1000;
    printf("%f ms.\n", tick);
}

void StopwatchPrintWithComment(const char* message, Stopwatch chrono) {
    printf(message, StopwatchGetElapsed(chrono));
}

void StopwatchPrintToFile(const char* text, Stopwatch chrono) {
    char buf[256];

    FILE *testfile = fopen("testfile.csv", "a");
    if (testfile == NULL) {
        printf("Impossible to load input files\n");
        return NULL;
    }
    snprintf(buf, sizeof buf, "%s,%s", text, StopwatchGetElapsed(chrono));
    fputs(buf, testfile);
    fclose(fopen);
}


#endif /* STOPWATCH_H_ */
