package sequential;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import model.Coordinate;

public class ClusterDNA {
	final private int MAX_SIM = 18;
	final private int DNA_LENGTH = 20;
	
	private int centroidNum;
	private String inputPath;
	private File inputFile;
	private ArrayList<char[]> strands;
	private ArrayList<char[]> centroids;
	private ArrayList<ArrayList<char[]>> clusters;

	public static void main(String[] args) {
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
		// step 4. Print clusters
		
	}

	public ClusterDNA(String inputPath, int centriodNum) {
		this.inputPath = inputPath;
		this.centroidNum = centriodNum;
		this.strands = new ArrayList<char[]>();
		this.centroids = new ArrayList<char[]>();
		this.clusters = new ArrayList<ArrayList<char[]>>();
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
			while ((line = br.readLine()) != null)
				strands.add(line.toCharArray());
			System.out.println("All sample strands are loaded.");
		} catch (IOException e) {
			throw new IOException();
		}
	}
	

	private void pickCentroids() throws Exception {
		int index = 0;
		centroids.add(strands.get(index++));
		while (centroids.size() < centroidNum) {
			if (index == strands.size()) 
				throw new Exception();
			if (isVaildCentroid(strands.get(index))) 
				centroids.add(strands.get(index));
			index++;
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
		while (true) {
			/* Classify all strands */
			for (char[] strand : strands) {
				int centroidIndex = 0;
				int maxSim = 0;
				for (int index = 0; index < centroidNum; index++) {
					if (centroids.get(index) == strand)
						continue;
					int curSim = sim(centroids.get(index), strand);
					if (curSim > maxSim) {
						maxSim = curSim;
						centroidIndex = index;
					}
				}
				clusters.get(centroidIndex).add(strand);
			}
			
			
			/* Calculate new centroids */
			
			int[][] countA = new int[centroidNum][DNA_LENGTH];
			int[][] countC = new int[centroidNum][DNA_LENGTH];
			int[][] countG = new int[centroidNum][DNA_LENGTH];
			int[][] countT = new int[centroidNum][DNA_LENGTH];
			
			/* Calculate sim for each strands in each cluster */
			for (int centroidIndex = 0; centroidIndex < centroidNum; centroidIndex++) {
				char[] centroid = centroids.get(centroidIndex);
				for (char[] strand : clusters.get(centroidIndex)) {
					for (int charIndex = 0; charIndex < DNA_LENGTH; charIndex++) {
						if (strand[charIndex] == centroid[charIndex]) {
							switch (strand[charIndex]) {
								case 'A':
									countA[centroidIndex][charIndex]++;
									break;
								case 'C':
									countC[centroidIndex][charIndex]++;
									break;
								case 'G':
									countG[centroidIndex][charIndex]++;
									break;
								case 'T':
									countT[centroidIndex][charIndex]++;
									break;
								default:
									break;
							}
						}
					}
				}
			}
			
			
			/* Find out new centroid for each cluster */
			for (int centroidIndex = 0; centroidIndex < centroidNum; centroidIndex++) {
				char[] newCentroid = new char[DNA_LENGTH];
				for (int charIndex = 0; charIndex < DNA_LENGTH; charIndex++) {
					
				}
			}
		}
	}
	
	private int sim(char[] strand1, char[] strand2) {
		int sim = 0;
		for (int i = 0; i < DNA_LENGTH; i++) {
			if (strand1[i] == strand2[i]) 
				sim++;
		}
		return sim;
	}

}
