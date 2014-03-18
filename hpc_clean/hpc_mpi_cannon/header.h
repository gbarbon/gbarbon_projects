/* 
 * File:   functions.h
 * Author: jian
 *
 * Created on 13 marzo 2014, 11.42
 */

#ifndef FUNCTIONS_H
#define	FUNCTIONS_H

//includes
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//define constants
#define TAG 13

//function declaration:
int* coordinate(int procNum, int totalProc);
void printmatrix();

#ifdef	__cplusplus
extern "C" {
#endif




#ifdef	__cplusplus
}
#endif

#endif	/* FUNCTIONS_H */

