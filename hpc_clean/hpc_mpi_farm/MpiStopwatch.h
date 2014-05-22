/*
 * Stopwatch.h
 *
 *  Created on: 08/mag/2013
 *      Author: henry
 */

#ifndef STOPWATCH_H_
#define STOPWATCH_H_

//#include <time.h>
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
    /*clock_t t1, t2;
    t1 = t2 = clock();

    // loop until t2 gets a different value
    while(t1 == t2)
            t2 = clock();*/

    double tick = MPI_Wtick() * 1000;
    printf("%f ms.\n", tick);
}

void StopwatchPrintWithComment(const char* message, Stopwatch chrono) {
    printf(message, StopwatchGetElapsed(chrono));
}

void StopwatchPrintToFile(const char* text, Stopwatch chrono) {
    char buf[256];
    char filename[256];

    /*locate home path*/
    char * homePath = getenv("HOME");
    snprintf(filename, sizeof filename, "%s/hpc_temp/hpc_time_res/testfile.csv", homePath);

    FILE *testfile = fopen(filename, "a");
    if (testfile == NULL) {
        printf("Impossible to load input files\n");
    } else {
        snprintf(buf, sizeof buf, "%s,%f\n", text, StopwatchGetElapsed(chrono));
        fputs(buf, testfile);
        fclose(testfile);
    }
}

void StopwatchPrintToFile2(const char* text, Stopwatch chrono) {
    char buf[256];
    char filename[256];

    /*locate home path*/
    char * homePath = getenv("HOME");
    snprintf(filename, sizeof filename, "%s/hpc_temp/hpc_time_res/testfile_farm.csv", homePath);

    FILE *testfile = fopen(filename, "a");
    if (testfile == NULL) {
        printf("Impossible to load input files\n");
    } else {
        snprintf(buf, sizeof buf, "%s,%f\n", text, StopwatchGetElapsed(chrono));
        fputs(buf, testfile);
        fclose(testfile);
    }
}

#endif /* STOPWATCH_H_ */
