import sys
import csv
import numpy
import getopt
import math
import random
import copy

def usage():
    print '$> python generatednadata.py <required args> [optional args]\n' + \
        '\t-c <#>\t\tNumber of clusters to generate\n' + \
        '\t-d <#>\t\tNumber of DNA strands per cluster\n' + \
        '\t-o <file>\tFilename for the output of the DNA data\n' + \
        '\t-l <#>\t\tLength of DNA strands\n'



def simDNA(strand1, strand2) :
	strandLen = len(strand1)
	sim = strandLen
	for i in range(strand1) :
		if strand1[i] != strand2[i] :
			sim = sim - 1
	return sim


def tooSimilar (strand, centroids, maxSim) :
	for dna in centroids :
		if sim(strand, dna) > maxSim :
			return False
	return True


def handleArgs(args):
	# set up return values
	numClusters = -1
	numDNA = -1
	output = None
	lenDNA = -1

	try:
		optlist, args = getopt.getopt(args[1:], 'c:p:v:o:')
	except getopt.GetoptError, err:
		print str(err)
		usage()
		sys.exit(2)

	for key, val in optlist:
		# first, the required arguments
		if   key == '-c':
			numClusters = int(val)
		elif key == '-d':
			numPoints = int(val)
		elif key == '-o':
			output = val
		# now, the optional argument
		elif key == '-l':
			lenDNA = int(val)

	# check required arguments were inputted  
	if numClusters < 0 or numPoints < 0 or maxValue < 1 or \
			output is None:
		usage()
		sys.exit()
	return (numClusters, numDNA, output, \
			lenDNA)




def drawOrigin(lenDNA):
	strand=[]
	tags = ['A','C','G','T']
	for i in range(lenDNA) :
		strand.append(tags[random.randint(0,3)])
	return strand 


# start by reading the command line
numClusters, numPoints, output, maxValue = handleArgs(sys.argv)

writer = csv.writer(open(output, "w"))

# step 1: generate each DNA centroid
centroids_radii = []
maxSim = 0.3 * lenDNA
for i in range(0, numClusters):
	centroid_radius = drawOrigin(lenDNA)
	# is it far enough from the others?
	while (tooSimilar(centroid_radius, centroids_radii, maxSim)):
		centroid_radius = drawOrigin(lenDNA)
	centroids_radii.append(centroid_radius)

# step 2: generate the points for each centroid
points = []
minClusterVar = 0.1 * lenDNA
maxClusterVar = 0.5 * lenDNA
for i in range(0, numClusters):
	# compute the variance for this cluster
	variance = numpy.random.uniform(minClusterVar, maxClusterVar)
	cluster = centroids_radii[i]
	for j in range(0, numDNA):
		acentroid=copy.deepcopy(cluster)
	# generate number of differences
	numDiff = int(abs(numpy.random.normal(0, deviation)))
	indDiff=numpy.random.randint(0, lenDNA, numDiff)
	for k in range(0,numDiff):
	    bases=['A','C','G','T']
	    origBase=acentroid[indDiff[k]]
	    bases.remove(origBase)
	    acentroid[indDiff[k]]=bases[random.randint(0,2)]
	# write the new strand out
	writer.writerow(acentroid)
	
	
