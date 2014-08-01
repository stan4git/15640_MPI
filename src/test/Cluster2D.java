package test;

import java.io.*;
import java.util.*;
import java.lang.Math;

/**
 * Cluster2D implements the k-means clustering algorithm. The program takes in a
 * list of 2D point locations, then separate them to target number of clusters.
 * 
 * @author Rui Zhang
 * @author Jing Gao
 * @version 1.0, 12/04/2013
 * @since 1.0
 */
public class Cluster2D {

	/**
	 * The main function of Cluster2D. It reads data file, chooses random
	 * initial centroids, does k-means until converge, and prints results.
	 * 
	 * @param argv
	 *            data file path and target number of clusters
	 * @since 1.0
	 */
	public static void main(String[] argv) throws Exception {

		if (argv.length != 2) {
			System.out
					.println("Proper Usage is: java Cluster2D datafile #cluster");
			System.exit(0);
		}
		if (!(new File(argv[0]).isFile())) {
			System.out.println("Please check input file path");
			System.exit(0);
		}
		try {
			int checkInt = Integer.parseInt(argv[1]);
			if (checkInt < 1) {
				System.out
						.println("Please provide correct number of target clusters");
				System.exit(0);
			}
		} catch (NumberFormatException nfe) {
			System.out
					.println("Please provide correct number of target clusters");
			System.exit(0);
		}

		/* the first argument is path of the data file containing the 2D points */
		String inputFile = argv[0];
		/* the second argument is the number of the clusters */
		int cenID = Integer.parseInt(argv[1]);

		/* scan input file and store the 2D points in an ArrayList */
		ArrayList<ArrayList<Float>> points = new ArrayList<ArrayList<Float>>();
		Scanner aScanner = new Scanner(new File(inputFile));
		String aLine = null;
		do {
			aLine = aScanner.nextLine();
			aLine = aLine.trim();
			if (aLine.length() > 1) {
				String[] parts = aLine.split(",");
				ArrayList<Float> xy = new ArrayList<Float>();
				xy.add(Float.parseFloat(parts[0]));
				xy.add(Float.parseFloat(parts[1]));
				points.add(xy);
			}
		} while (aScanner.hasNext());

		/* choose initial centroids, use either one of the two strategies below */

		/*
		 * strategy 1 randomly select initial centroids from the original 2D
		 * points, and make sure the centroids are far away from each other
		 */
		float minDist = 0.5f; // set it smaller than minDistance in
								// generaterawdata.py
		ArrayList<ArrayList<Float>> centroids = new ArrayList<ArrayList<Float>>();
		int randInd = (int) (Math.random() * points.size());
		centroids.add(points.get(randInd));
		for (int i = 0; i < cenID; i++) {
			while (tooClose(points.get(randInd), centroids, minDist)) {
				randInd = (int) (Math.random() * points.size());
			}
			centroids.add(points.get(randInd));
		}

		/*
		 * strategy 2 use the first k points as initial centroids
		 */
		/*
		 * ArrayList<ArrayList<Float>> centroids = new
		 * ArrayList<ArrayList<Float>>(); for (int i=0;i<k_param;i++){
		 * centroids.add(points.get(i)); }
		 */

		/* the threshold to determine if centroids are stable */
		double threshold = 1.0e-6;
		/* the total difference of centroids from two iterations */
		double incre = 1.0e+6;
		/* iterate until centroids are stable */
		while (incre > threshold) {

			/* record sum of x & y for each cluster */
			float[][] tagXY = new float[cenID][2];
			/* record size of each cluster */
			int[] tagSize = new int[cenID];
			/* record tag of each point */
			ArrayList<Integer> tags = new ArrayList<Integer>();

			/*
			 * assign tag for each point according to minimum distance among
			 * centroids
			 */
			for (int i = 0; i < points.size(); i++) {
				ArrayList<Float> dists = new ArrayList<Float>();
				for (int j = 0; j < cenID; j++) {
					dists.add(dist(points.get(i), centroids.get(j)));
				}
				int tag = dists.indexOf(Collections.min(dists));
				tags.add(tag);
				tagXY[tag][0] += points.get(i).get(0);
				tagXY[tag][1] += points.get(i).get(1);
				tagSize[tag] += 1;
			}

			/*
			 * compute the total difference of centroids from last iteration,
			 * see if converged
			 */
			incre = 0.0;
			for (int i = 0; i < cenID; i++) {
				tagXY[i][0] /= tagSize[i];
				tagXY[i][1] /= tagSize[i];
				double tempX2 = Math.pow(
						(double) (tagXY[i][0] - centroids.get(i).get(0)), 2);
				double tempY2 = Math.pow(
						(double) (tagXY[i][1] - centroids.get(i).get(1)), 2);
				incre = incre + Math.sqrt(tempX2 + tempY2);
				ArrayList<Float> updatedCentroid = new ArrayList<Float>();
				updatedCentroid.add(tagXY[i][0]);
				updatedCentroid.add(tagXY[i][1]);
				centroids.set(i, updatedCentroid);
			}

			/* if converged, print out the results */
			if (incre <= threshold) {
				for (int i = 0; i < cenID; i++) {
					System.out.println("centroid " + i + ":");
					System.out.println(centroids.get(i));
					System.out.println("points:");
					for (int j = 0; j < points.size(); j++) {
						if (tags.get(j) == i) {
							System.out.println(points.get(j));
						}
					}
				}
			}

			System.out.println(incre);
		}

	}

	/*
	 * compute Euclidean distance between two points
	 * 
	 * @param p1 the 2D point
	 * 
	 * @param p2 the other 2D point
	 * 
	 * @return the Euclidean distance between the two points
	 * 
	 * @since 1.0
	 */
	public static Float dist(ArrayList<Float> p1, ArrayList<Float> p2) {
		double xDiff = p1.get(0) - p2.get(0);
		double yDiff = p1.get(1) - p2.get(1);
		float dis = (float) Math.sqrt(Math.pow(xDiff, 2) + Math.pow(yDiff, 2));
		return dis;
	}

	/*
	 * decide if one point is too close to any of the points in a list
	 * 
	 * @param one the 2D point
	 * 
	 * @param list a list of 2D points to be compared with
	 * 
	 * @param minDist the min Euclidean distance allowed
	 * 
	 * @return true if the point is too close to one of the points in list
	 * 
	 * @since 1.0
	 */
	public static boolean tooClose(ArrayList<Float> one,
			ArrayList<ArrayList<Float>> list, Float minDist) {
		for (int iter = 0; iter < list.size(); iter++) {
			if (dist(one, list.get(iter)) < minDist) {
				return true;
			}
		}
		return false;
	}

}
