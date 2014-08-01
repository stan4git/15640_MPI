package sequential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;

public class ClusterDNA {
	final private int MAX_SIM = 18;
	final private int DNA_LENGTH = 20;
	
	private int centroidNum;
	private String inputPath;
	private File inputFile;
	private ArrayList<char[]> strands;
	private ArrayList<char[]> centroids;
//	private ArrayList<ArrayList<char[]>> clusters;

	public static void main(String[] args) {
		long startTime = System.currentTimeMillis();
		/* java ClusterDNA <input> <clusterNum> */
		// step 1. Load arguments
		if (args.length != 2) {
			System.err.println("Wrong number of arguments. "
					+ "Usage: \"java cluster2D <input_file> <centriod_num>\".");
			System.exit(-1);
		}
		String inputPath = null;
		int centriodNum = 0;
		ClusterDNA cluster = null;
		try {
			inputPath = args[0];
			centriodNum = Integer.parseInt(args[1]);
			cluster = new ClusterDNA(inputPath, centriodNum);
			cluster.loadDNA();
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

		// step 3. Find clusters
		cluster.findCluster();
		
		// step 4. Print clusters
		cluster.printCluster();
		
		System.out.println("Running time: " + ((System.currentTimeMillis() - startTime) / 1000.00) + " seconds.");
	}

	public ClusterDNA(String inputPath, int centriodNum) {
		this.inputPath = inputPath;
		this.centroidNum = centriodNum;
		this.strands = new ArrayList<char[]>();
		this.centroids = new ArrayList<char[]>();
	}
	

	private void loadDNA() throws IOException {
		this.inputFile = new File(inputPath);
		if (!inputFile.exists()) {
			System.err.println("File not found.");
			throw new FileNotFoundException();
		}
		BufferedReader br = new BufferedReader(new FileReader(inputFile));
		String line;
		try {
			while ((line = br.readLine()) != null) {
				line = line.trim().replace(",", "");
				if (line.length() > 1) 
					strands.add(line.toCharArray());
			}
			System.out.println("All sample points are loaded.");
		} catch (IOException e) {
			throw new IOException();
		} finally {
			br.close();
		}
	}
	

	private void pickCentroids() throws Exception {
		/* Validate distances between centroids. */
		int index = 0;
		while (centroids.size() < centroidNum) {
			do {
				index = getRamdonPointIndex(strands.size());
			} while (!isVaildCentroid(strands.get(index)));
			centroids.add(strands.get(index));
		}
	}
	
	
	private boolean isVaildCentroid(char[] strand) {
		for (char[] centroid : centroids) {
			if (sim(strand, centroid) > MAX_SIM) 
				return false;
		}
		return true;
	}
	
	
	private void findCluster() {
		boolean converged = false;
		char[] tagList = new char[] {'A','C','T','G'};
		int count = 1;
		while (!converged) {
			int[][] countA = new int[centroidNum][DNA_LENGTH];
			int[][] countC = new int[centroidNum][DNA_LENGTH];
			int[][] countG = new int[centroidNum][DNA_LENGTH];
			int[][] countT = new int[centroidNum][DNA_LENGTH];
			
			/* Classify all strands */
			for (char[] strand : strands) {
				ArrayList<Integer> simList = new ArrayList<Integer>(centroidNum);
				/* Calculate sim for each strands in each cluster */
				for (int cenIndex = 0; cenIndex < centroidNum; cenIndex++) {
					int sim = sim(centroids.get(cenIndex), strand);
					simList.add(sim);
				}
				/* Find the most similar centroid */
				int selectedCentroid = simList.indexOf(Collections.max(simList));
				
				/* Add apparance to record */
				for (int charIndex = 0; charIndex < DNA_LENGTH; charIndex++) {
					switch (strand[charIndex]) {
					case 'A':
						countA[selectedCentroid][charIndex]++;
						break;
					case 'C':
						countC[selectedCentroid][charIndex]++;
						break;
					case 'G':
						countG[selectedCentroid][charIndex]++;
						break;
					case 'T':
						countT[selectedCentroid][charIndex]++;
						break;
					default:
						break;
					}
				}
			}
			
			
			
			/* Find out new centroid for each cluster */
			ArrayList<char[]> newCentroids = new ArrayList<char[]>();
			
			for (int centroidIndex = 0; centroidIndex < centroidNum; centroidIndex++) {
				char[] centroid = new char[DNA_LENGTH];
				for (int charIndex = 0; charIndex < DNA_LENGTH; charIndex++) {
					ArrayList<Integer> simArr = new ArrayList<Integer>();
					simArr.add(countA[centroidIndex][charIndex]);
					simArr.add(countC[centroidIndex][charIndex]);
					simArr.add(countT[centroidIndex][charIndex]);
					simArr.add(countG[centroidIndex][charIndex]);
					centroid[charIndex] = tagList[simArr.indexOf(Collections.max(simArr))];
				}
				newCentroids.add(centroid);
			}
			
			
			/* If converged, replace the old centroids */
			if (isConverged(newCentroids)) {
				converged = true;
			} else {
				System.out.println(count++);
			}
			centroids = newCentroids;
		}
	}
	
	
	private boolean isConverged(ArrayList<char[]> newCentroids) {
		if (newCentroids == null || newCentroids.size() != centroidNum) 
			return false;
		for (int centroidIndex = 0; centroidIndex < centroidNum; centroidIndex++) 
			for (int charIndex = 0; charIndex < DNA_LENGTH; charIndex++) 
				if (newCentroids.get(centroidIndex)[charIndex] != centroids.get(centroidIndex)[charIndex]) 
					return false;
		return true;
	}
	
	private int getRamdonPointIndex(int max) {
		return (int)(Math.random() * max);
	}
	
	
	private int sim(char[] strand1, char[] strand2) {
		int sim = 0;
		for (int i = 0; i < DNA_LENGTH; i++) {
			if (strand1[i] == strand2[i]) 
				sim++;
		}
		return sim;
	}
	
	private void printCluster() {
		for (int centroidIndex = 0; centroidIndex < centroidNum; centroidIndex++) {
			System.out.println("Cluster " + centroidIndex + ":");
			String str = new String(centroids.get(centroidIndex));
			System.out.println("Centroid: " + str);
		}
	}
}
