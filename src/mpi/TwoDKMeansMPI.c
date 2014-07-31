#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include "Util.h"

/* This is the threshold to end the calculation for 2D */
#define TwoD_DIFF_THRESHOLD 0.0001
/* This is the max difference of 2D */
#define MAX_DIFF 2147483647

int main(int argc,char** argv){
    
    /* Step 1: parsing the user's input */
    
    /* total lines of content */
    int lineNums;
    
    /* file name */
    char* filename;
    
    /* the cluster numbers */
    int cluster;
    
    /* dimension of each point or DNA squence */
    int dimension;
    
    /* check the arguments */
	if (argc < 4) {
		printf("Usage: TwoDKMeansMPI <input file name> <line Numbers> <cluster Numbers> \n");
		exit(-1);
	}
    
    /* file name */
	filename = argv[1];
    
	/* how many points */
	lineNums = atoi(argv[2]);
    
	/* how many clusters */
	cluster = atoi(argv[3]);
    
    /* 2D dimension */
    dimension = 2;
    
    
    
    /* Step 2: Initilize the MPI */
    
    /* Process number, Process ID */
    int numprocs, rank;
    
    /* Initilize of MPI*/
    MPI_Init(&argc,&argv);
    /*declare a variable to hold the time returned*/
    double startwtime, endwtime;
    /*get the time just before work to be timed*/
    startwtime = MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    
    
    
    
    /* Step 3: Build the buffers */
    
    /* the total contents to be handled for 2D points */
    double* TwoDSource;
    
    /* the contents for each processor of 2D points */
    double* Recvbuf2D;
    
    /* 2D centroids */
    double* TwoDCentroids;
    
    /* the working point numbers of each process */
    
    int handleNumbers;
    
    /* for example: 501 lines and 10 processor
     * it need 51 lines per processor, but if it just has 500 lines,
     * it need 50 lines per processor.
     */
    if (lineNums % numprocs > 0) {
        handleNumbers = (lineNums / numprocs  + 1 ) * dimension;
    } else {
        handleNumbers = lineNums / numprocs * dimension;
    }
    
    TwoDSource = malloc(sizeof(double) * dimension * lineNums);
    Recvbuf2D = malloc(sizeof(double) * handleNumbers);
    TwoDCentroids = malloc(sizeof(double) * cluster * dimension);
    
    
    
    
    /* Step 4: Read the contents, generating centroids and calculating the cursor */
    
    /* the file to be handled */
    FILE* fp;
    
    /* if it's master, read the data source and genterate centroids */
	if (rank == 0) {
		fp = fopen(filename, "r");
		if(fp == NULL) {
			printf("Cannot open file %s\n", filename);
			exit(-1);
		}
		
		read2DContents(fp, TwoDSource);
        generate2DCentroids(TwoDCentroids, TwoDSource, lineNums, dimension, cluster);
        
		fclose(fp);
	}
    
    /* compute the handling lines and start index for each processes */
	int* sendcounts = malloc(sizeof(int)*numprocs);
	int* displs = malloc(sizeof(int)*numprocs);
	int processIndex;
    /* the previous n - 1 process must be full */
	for(processIndex = 0;processIndex < numprocs - 1;processIndex++) {
		sendcounts[processIndex] = handleNumbers;
	}
    /* tha last process's handling lines */
	sendcounts[numprocs - 1] = lineNums * dimension - handleNumbers * (numprocs - 1);
    
    
    /* Compute the beginning index of line number for each process */
	int offset = 0;
    /* the two variable are used in the loop */
    int i,j;
	for(i = 0;i < numprocs;i++) {
		displs[i] = offset;
		offset += sendcounts[i];
	}
    
	/* scatter the data and broadcast the centroids*/
	MPI_Scatterv (TwoDSource,sendcounts,displs,MPI_DOUBLE,Recvbuf2D,handleNumbers,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast (TwoDCentroids,cluster * dimension,MPI_DOUBLE,0,MPI_COMM_WORLD);
    
    
    
    
    /* Step5: K-Means Calculating */
    
    /* categorized all the points or DNA strands into different clusters for each processor*/
    int handleRows = sendcounts[rank] / dimension;
    int* categories = malloc(sizeof(int) * handleRows);
    /* all the labels of all the points on the master process */
	int* totalCategories = malloc(sizeof(int)*lineNums);
    
    do {
        /* the flag for temination of the while loop */
        int flag = 0;
        /* the distributed centroids from master node */
        double* distributedCentroids = malloc(sizeof(double) * cluster * dimension);
        /* the new generated centroids */
        double* newGeneratedCentroids = malloc(sizeof(double) * cluster * dimension);
        memset(newGeneratedCentroids, 0, sizeof(double) * cluster * dimension);
        /* the numbers of points in distributed contents */
        int* distributedPoints = malloc(sizeof(int) * cluster);
        /* the numbers of new generated points in each processor */
        int* newGeneratedPoints = malloc(sizeof(int) * cluster);
        memset(newGeneratedPoints, 0, sizeof(int) * cluster);
        
        /* Calculate the category and new Centroids' sum */
        for(i = 0; i < handleRows; i++) {
            int category = -1;
            double TwoDmaxDiff = MAX_DIFF;
            for(j = 0; j < cluster; j++) {
                double calculatedDistance = TwoDDistance(TwoDCentroids + j * dimension, Recvbuf2D + i * dimension);
                if (calculatedDistance < TwoDmaxDiff) {
                    TwoDmaxDiff =calculatedDistance;
                    category = j;
                }
                categories[i] = category;
                newGeneratedPoints[category]++;
                /* sum each point */
                for (j = 0; j < dimension; j++) {
                    newGeneratedCentroids[category * dimension + j] += Recvbuf2D[i * dimension + j];
                }
            }
        }
        /* reduce step - new centroid sum to the master process */
        MPI_Reduce(newGeneratedCentroids, distributedCentroids, cluster * dimension, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        /* reduce step - the cluster counts to the master process */
        MPI_Reduce(newGeneratedPoints, distributedPoints, cluster, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        
        /* handle the results from different processor */
        if (rank == 0){
            /* intilize the difference sum */
            double sumDistance = 0;
            
            /* calculate new centroids' average value */
            for (i = 0; i < cluster; i++) {
                for (j = 0; j < dimension; j++) {
                    distributedCentroids[i * dimension + j] /= distributedPoints[i];
                }
            }
            
            /* iterate each cluster and update the sum */
            for (i = 0; i < cluster; i++) {
                sumDistance += TwoDDistance(distributedCentroids + i * dimension,TwoDCentroids + i * dimension);
            }
            /* if the difference is less than a threshold, set the termination flag */
            if (sumDistance < TwoD_DIFF_THRESHOLD) {
                flag = 1;
            }
            /* free the TwoDCentroids and give the TwoDCentroids new values */
            free(TwoDCentroids);
            TwoDCentroids = distributedCentroids;
        }
        /* broadcast the temination flag */
        MPI_Bcast (&flag,1,MPI_INT,0,MPI_COMM_WORLD);
        
        free(newGeneratedPoints);
        free(distributedPoints);
        free(newGeneratedCentroids);
        
        if(flag) {
            break;
        }
        
        /* broadcast the new centroids */
        MPI_Bcast (TwoDCentroids,cluster * dimension,MPI_DOUBLE,0,MPI_COMM_WORLD);
    }while (1);
    
    /* gather all the labels of all the points on the master process */
	int displacement = 0;
	for(i = 0;i < numprocs;i++) {
		sendcounts[i] /= dimension;
		displs[i] = displacement;
		displacement += sendcounts[i];
	}
	MPI_Gatherv (categories,handleRows,MPI_INT,totalCategories,sendcounts,displs,MPI_INT,0,MPI_COMM_WORLD);
	
	
    
    
    /* Step6: GC */
	free(sendcounts);
	free(displs);
	free(categories);
	free(totalCategories);
    free(TwoDSource);
    free(Recvbuf2D);
    free(TwoDCentroids);
    
    /*get the time just after work is done and take the difference */
    endwtime = MPI_Wtime();
    printf("Timing span of this job is %lf seconds.\n",endwtime-startwtime);
	MPI_Finalize();
}
