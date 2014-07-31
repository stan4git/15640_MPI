package sequential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

public class Cluster2D {
	final double FLUC_RANGE = 1.0e-3;
	final double MIN_CEN_DIST = 0.3;
	
	private int centroidNum;
	private String inputPath;
	private File inputFile;
	private ArrayList<Coordinate> points;
	private ArrayList<Coordinate> centroids;
	private ArrayList<ArrayList<Coordinate>> clusters;
//	private ArrayList<Integer> assignCluster;
	
	
	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		// "java cluster2D <input_file> <centriod_num>"
		// step 1. Load arguments
		if (args.length != 2) {
			System.err.println("Wrong number of arguments. "
					+ "Usage: \"java cluster2D <input_file> <centriod_num>\".");
			System.exit(-1);
		}
		String inputPath = null;
		int centriodNum = 0;
		Cluster2D cluster = null;
		try {
			inputPath = args[0];
			centriodNum = Integer.parseInt(args[1]);
			cluster = new Cluster2D(inputPath, centriodNum);
			
			cluster.loadPoints();
		// step 2. Pick initial centroids
			cluster.pickCentroids();
		} catch (NumberFormatException e1) {
			System.err.println("Invaild centroid number. Please try again.");
			System.exit(-1);
		} catch (IOException e2) {
			System.err.println("IOException occurs when initilizing. Please try again.");
			System.exit(-1);
		} catch (Exception e3) {
			System.err.println("Cannot find enough centroids. Please validate input data.");
			System.exit(-1);
		}

		// step 3. Calculate distance and assign points to centroids
		// step 4. Calculate new centroids
		// step 5. Check fluctuation, loop if fail
		cluster.findClusters();
		
		// step 6. Output centroids
		cluster.printClusters();
		System.out.println("Running time: " + ((System.currentTimeMillis() - startTime) / 1000.00) + " seconds.");
	}

	
	
	
	
	public Cluster2D(String inputPath, int centroidNum) {
		this.centroidNum = centroidNum;
		this.inputPath = inputPath;
		this.points = new ArrayList<Coordinate>();
		this.centroids = new ArrayList<Coordinate>();
		initCluster();
	}

	
	
	
	
	/*
	 * This method is used to load all the sample points into
	 * memory.
	 */
	public void loadPoints() throws IOException {
		this.inputFile = new File(inputPath);
		if (!inputFile.exists()) {
			System.err.println("File not found.");
			throw new FileNotFoundException();
		}
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		String line;
		try {
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.length() > 1) {
					String[] coor = line.split(",");
					float x = Float.parseFloat(coor[0]);
					float y = Float.parseFloat(coor[1]);
					Coordinate point = new Coordinate(x, y);
					points.add(point);
				}
			}
			System.out.println("All sample points are loaded.");
		} catch (IOException e) {
			throw new IOException();
		} finally {
			br.close();
		}
	}
	
	
	
	
	
	/*
	 * Used to pick out designated number of centroids.
	 * To ensure these centroids are not too close, isTooClose will be invoked
	 * and MIN_CEN_DIST will be used to define the minimum distance.
	 */
	private void pickCentroids() throws Exception {
		//Pick first point as centroids.
		int cursor = 0;
//		centroids.add(points.get(cursor++));
		//Validate distances between centroids.
		while (centroids.size() < centroidNum) {
			if (cursor == points.size()) 
				throw new Exception();
			if (cursor == 0) {
				centroids.add(points.get(getRamdonPointIndex(points.size())));
			} else {
				if (isValidCentroid(points.get(cursor))) 
					centroids.add(points.get(cursor));
			}
			cursor++;
		}
	}
	
	
	private boolean isValidCentroid(Coordinate currentPoint) {
		for (Coordinate centroid : centroids) {
			if (dist(currentPoint, centroid) < MIN_CEN_DIST)
				return false;
		}
		return true;
	}

	
	private int getRamdonPointIndex(int max) {
		return (int)Math.random() * max;
	}
	
	
	private void initCluster() {
		this.clusters = new ArrayList<ArrayList<Coordinate>>();
		for (int i = 0; i < centroidNum; i++) {
			clusters.add(new ArrayList<Coordinate>());
		}
//		assignCluster = new ArrayList<Integer>();
	}
	
	private void findClusters() {
		float fluc = Float.MAX_VALUE;
		while (fluc > FLUC_RANGE) {
			//generate a new cluster list
			initCluster();
			for (int pointNum = 0; pointNum < points.size(); pointNum++) {
//				if (centroids.contains(point))
//					continue;
				int clusterIndex = 0;
				double minDist = dist(points.get(pointNum), centroids.get(0));
				for (int i = 1; i < centroidNum; i++) {
					Coordinate curCentroid = centroids.get(i);
					double curDist = dist(points.get(pointNum), curCentroid);
					if (curDist < minDist) {
						clusterIndex = i;
						minDist = curDist;
					}
				}
				clusters.get(clusterIndex).add(points.get(pointNum));
//				assignCluster.set(pointNum, clusterIndex);
			}
			
			//calculate new centroids
			ArrayList<Coordinate> newCentroids = new ArrayList<Coordinate>();
			for (ArrayList<Coordinate> cluster : clusters) {
				newCentroids.add(calculateMedian(cluster));
			}
			//calculate fluctuation
			fluc = calculateFluctuation(newCentroids);
			centroids = newCentroids;
			System.out.println(fluc);
		}
	}
	
	
	private Coordinate calculateMedian(ArrayList<Coordinate> cluster) {
		int size = cluster.size();
		float xSum = 0;
		float ySum = 0;
		for (Coordinate point : cluster) {
			xSum += point.x;
			ySum += point.y;
		}
		Coordinate newCentroid = new Coordinate(xSum / size, ySum / size);
		return newCentroid;
	}
	
	
	private float calculateFluctuation(ArrayList<Coordinate> newCentroids) {
		float fluc = 0;
		for (int i = 0; i < centroidNum; i++) {
			fluc += dist(newCentroids.get(i), centroids.get(i));
		}
		return fluc;
	}
	
	
	private void printClusters() {
		for (int index = 0; index < centroidNum; index++) {
			System.out.println("Cluster " + index + ": ");
			System.out.println("Centroid: " + centroids.get(index).toString());
//			System.out.print("Points: ");
//			for (Coordinate point : clusters.get(index)) {
//				System.out.print(point.toString() + " ");
//			}
			System.out.println();
		}
	}
	
	
	private double dist(Coordinate point1, Coordinate point2) {
		double dist = 0d;
		dist = Math.sqrt(Math.pow(point1.x - point2.x, 2) + Math.pow(point1.y - point2.y, 2));
		return dist;
	}
}
