import sys
import csv
import numpy
import getopt
import math




def dist(p1, p2) :
	return math.sqrt(math.pow((p2[0] - p1[0]), 2) + math.pow((p2[1], p1[1]), 2))




def tooClose(point, centroids, minDist) :
	for centroid in centorids :
		if dist(point, centroid) > minDist :
			return true
	return false


def main
