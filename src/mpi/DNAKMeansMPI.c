#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include "Util.h"

/* This is the threshold to end the calculation for DNA */
#define DNA_DIFF_THRESHOLD 1
/* This is the max difference of DNA */
#define MAX_DIFF 2147483647
/* This is the basic components fo the DNA sequence */
char DNA_Components[4] = {'A','C','G','T'};

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
    
	/* check the arguments if it is the DNA case */
	if (argc < 5) {
		printf("Usage: DNAKMeansMPI <input file name> <line Numbers> <cluster Numbers> <dimension>\n");
		exit(-1);
	}
    
    /* file name */
	filename = argv[1];
    
	/* how many points */
	lineNums = atoi(argv[2]);
    
	/* how many clusters */
	cluster = atoi(argv[3]);
    
    /* decide the dimension */
    dimension = atoi(argv[4]);
    

    
    
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
    
    /* the total contents to be handled for DNA strands */
    char* DNASource;
    
    /* the contents for each processor of DNA strands */
    char* RecvbufDNA;
    
    /* DNA centroids */
    char* DNACentroids;
    
    /* the working line numbers of each process */
    
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
    
    DNASource = malloc(sizeof(char) * dimension * lineNums);
    RecvbufDNA = malloc(sizeof(char) * handleNumbers);
    DNACentroids = malloc(sizeof(char) * cluster * dimension);
    
    
    
    
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
		
		readDNAContents(fp, DNASource, dimension);
        generateDNACentroids(DNACentroids, DNASource, lineNums, dimension, cluster);
        
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
    
	/* scatter the data and broadcast the center*/
	MPI_Scatterv (DNASource,sendcounts,displs,MPI_CHAR,RecvbufDNA,handleNumbers,MPI_CHAR,0,MPI_COMM_WORLD);
    MPI_Bcast (DNACentroids,cluster * dimension,MPI_CHAR,0,MPI_COMM_WORLD);
    
    
    
    
    /* Step5: K-Means Calculating */
    
    /* categorized all the points or DNA strands into different clusters for each processor*/
    int handleRows = sendcounts[rank] / dimension;
    int* categories = malloc(sizeof(int) * handleRows);
    /* all the labels of all the points on the master process */
	int* totalCategories = malloc(sizeof(int)*lineNums);
    
    do {
            /* termination flag */
			int flag = 0;
            
            /* Example: Cluster = 2, Dimension = 4 */
            /* A             C              G               T          */
            /*  0| 4| 8|12   1 | 5| 9|13    2| 6|10|14      3| 7|11|15 */
            /* 16|20|24|28   17|21|25|29   18|22|26|30     19|23|27|31 */
             
			int* newGeneratedContents = malloc(sizeof(int) * cluster * dimension * 4);
			memset(newGeneratedContents, 0, sizeof(int) * cluster * dimension * 4);
			int* distributedContents = malloc(sizeof(int) * cluster * dimension * 4);
            
			/* the centroids for every cluster computed */
			char* distributedCentroids = malloc(sizeof(char) * cluster * dimension);
            
			/* Calculate each points distances to centroids and categorize it */
			for(i = 0; i < handleRows; i++) {
				int category = -1;
				int DNAmaxDiff = MAX_DIFF;
				for(j = 0; j < cluster;j++) {
					int calculatedDistance = DNADistance(DNACentroids+j*dimension, RecvbufDNA+i*dimension, dimension);
					if (calculatedDistance < DNAmaxDiff) {
						DNAmaxDiff = calculatedDistance;
						category = j;
					}
				}
                
				/* Work for each Line's character and update the array's value */
				for(j = 0; j < dimension; j++) {
					switch(RecvbufDNA[i * dimension + j]) {
						case 'A':
							newGeneratedContents[4 * dimension * category + 4 * j + 0]++;
							break;
						case 'C':
							newGeneratedContents[4 * dimension * category + 4 * j + 1]++;
							break;
						case 'G':
							newGeneratedContents[4 * dimension * category + 4 * j + 2]++;
							break;
						case 'T':
							newGeneratedContents[4 * dimension * category + 4 * j + 3]++;
							break;
						default:
							break;
					}
				}
                
				/* update the point's category */
				categories[i] = category;
			}
			
			/* reduce step - the character counts in every position of each point on the master node */
			MPI_Reduce(newGeneratedContents, distributedContents, cluster * 4 * dimension, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            
            
			
			/* on the master node, compute the new centroid for each cluster */
			if(rank == 0) {
                /* i means cluster */
				for (i = 0; i < cluster; i++) {
					/* j means each character of one sequence */
					for (j = 0; j < dimension; j++) {
						int aValue = distributedContents[i * dimension * 4 + 4 * j + 0];
                        int cValue = distributedContents[i * dimension * 4 + 4 * j + 1];
                        int gValue = distributedContents[i * dimension * 4 + 4 * j + 2];
                        int tValue = distributedContents[i * dimension * 4 + 4 * j + 3];
                        int resultValue = max(aValue,max(cValue,max(gValue,tValue)));
                        if (resultValue == aValue) {
                            distributedCentroids[i * dimension + j] = DNA_Components[0];
                        } else if(resultValue == cValue) {
                            distributedCentroids[i * dimension + j] = DNA_Components[1];
                        } else if (resultValue == gValue) {
                            distributedCentroids[i * dimension + j] = DNA_Components[2];
                        } else {
                            distributedCentroids[i * dimension + j] = DNA_Components[3];
                        }
					}
				}
				/* initialize the diffrenciation of centroids between orginal and new genereated */
				int sumDistance = 0;
				/* iterate each cluster centroid and update the difference of centroids */
				for(i =0;i<cluster;i++) {
					sumDistance += DNADistance(distributedCentroids + dimension * i,DNACentroids + dimension * i,dimension);
				}
				/* see if it is time to terminate */
				if(sumDistance <= DNA_DIFF_THRESHOLD) {
                    flag = 1;
                }
				free(DNACentroids);
				DNACentroids = distributedCentroids;
			}
            
			/* broadcast the terminate flag */
			MPI_Bcast (&flag,1,MPI_INT,0,MPI_COMM_WORLD);
			free(newGeneratedContents);
			free(distributedContents);
			if(flag) {
                break;
            }
			
            /* broadcast the new centroids */
			MPI_Bcast (DNACentroids,cluster * dimension,MPI_CHAR,0,MPI_COMM_WORLD);
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
    free(DNASource);
    free(RecvbufDNA);
    free(DNACentroids);
    
    /*get the time just after work is done and take the difference */
    endwtime = MPI_Wtime();
    printf("Timing span of this job is %lf seconds.\n",endwtime - startwtime);
	MPI_Finalize();
}
