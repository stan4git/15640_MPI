#ifndef _UTIL_H_
#define _UTIL_H_

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "string.h"

/*
 * read 2D file
 */
void read2DContents(FILE* file, double* array);
/*
 * read DNA file
 */
void readDNAContents(FILE* file, char* array, int dimension);

/*
 * generate 2D centroids
 */
void generate2DCentroids(double* centroids, double* source, int lineNums, int dimension, int cluster);

/*
 * generate DNA centroids
 */
void generateDNACentroids(char* centroids, char* source, int lineNums, int dimension, int cluster);

/*
 * 2D : compute the distance
 */
double TwoDDistance(double* centroid, double* point);

/*
 * DNA : compute the similarity two DNA strands
 */
int DNADistance(char* centroid, char* strand, int dimension);

/*
 * determin whether the source point are too close to one of the centroids
 */
int tooClose(double* centroids, double* source, int num);

/*
 * determin whether the source strand are too similar to one of the centroids
 */
int tooSimilar(char* centroids, char* source, int dimension, int num);

/*
 * max value between two values
 */
int max(int num1,int num2);

#endif

