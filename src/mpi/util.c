#include "Util.h"

/*
 * read 2D file
 */
void read2DContents(FILE* file, double* array) {
	int index = 0;
	/* file not open */
	if (file == NULL) {
		printf("File not open\n");
		exit(1);
	}
	
	/* scan the whole file line by line */
	while(fscanf(file,"%lf,%lf\n",array + 2 * index,array + 2 * index + 1) != EOF) {
		index++;
	}
	
}

/*
 * read DNA file
 */
void readDNAContents(FILE* file, char* array, int dimension) {
    
	if (file == NULL) {
		printf ("File not open\n");
		exit(1);
	}
	
	int index = 1;
	while (1) {
		int result;
		if(index % dimension == 0) {
			result = fscanf(file,"%c\n", array + index - 1);
		} else {
			result = fscanf(file,"%c,", array + index - 1);
		}
		if (result == EOF)
			break;
		index++;
	}
}


/*
 * generate 2D centroids
 */
void generate2DCentroids(double* centroids, double* source, int lineNums, int dimension, int cluster) {
	int i;
    /* select the random line */
	int index = ((int)rand()) % lineNums;
    /* copy the select line into the centroids array */
	memcpy(centroids, source + index * dimension,dimension * sizeof(double));
	for (i = 0;i < cluster;i++) {
		/* if the selected point is too close with any selected centroids,
         * it will reselect the point */
        while(tooClose(centroids, source + index * dimension, i)) {
			index = ((int)rand())%lineNums;
		}
		memcpy(centroids + i*dimension, source+index*dimension,dimension * sizeof(double));
	}
}

/*
 * generate DNA centroids
 */
void generateDNACentroids(char* centroids, char* source, int lineNums, int dimension, int cluster) {
	int i;
    /* select the random line */
	int index = ((int)rand())%lineNums;
    /* copy the select line into the centroids array */
	memcpy(centroids,source+index*dimension,dimension);
	index = ((int)rand() )% lineNums;
	for (i = 1;i < cluster;i++) {
		while(tooSimilar(centroids, source+index*dimension, dimension, i)) {
			index = ((int)rand())%lineNums;
		}
		memcpy(centroids+i*dimension, source+index*dimension, dimension);
	}
}


/*
 * 2D : compute the distance
 */
double TwoDDistance(double* centroid, double* point) {
    return sqrt(pow(centroid[0]-point[0],2) + pow(centroid[1] - point[1],2));
}

/*
 * DNA : compute the similarity two DNA strands
 */
int DNADistance(char* centroid, char* strand, int dimension) {
    int i;
    int sum = dimension;
    for(i = 0;i<dimension;i++){
        if (centroid[i] == strand[i]) {
            sum--;
        }
    }
    return sum;
}

/*
 * determin whether the source point are too close to one of the centroids
 */
int tooClose(double* centroids, double* source, int num) {
	double minDist = 0.5;
	int i;
	for(i = 0;i<num;i++){
		if(TwoDDistance(centroids+i*2, source) < minDist)
			return 1;
	}
	return 0;
}

/*
 * determin whether the source strand are too similar to one of the centroids
 */
int tooSimilar(char* centroids, char* source, int dimension, int num) {
	int minDist = (int)(0.3 * dimension);
	int i;
	for(i = 0;i<num;i++){
		if(DNADistance(centroids+i*dimension, source, dimension) < minDist)
			return 1;
	}
	return 0;
	
}


/*
 * max value between two values
 */
int max(int num1,int num2) {
    int result;
    
    if (num1 > num2)
        result = num1;
    else
        result = num2;
    
    return result;
}



