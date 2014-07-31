package test;

import java.io.*;
import java.util.*;
import java.lang.Math;

/**
 * ClusterDNA implements the k-means clustering algorithm.
 * The program takes in a list of DNA strands,
 * then separate them to target number of clusters.
 *
 * @author      Rui Zhang
 * @author      Jing Gao
 * @version     1.0, 12/04/2013
 * @since       1.0
 */
public class ClusterDNA {

    /**
     * The main function of ClusterDNA.
     * It reads data file, chooses random initial centroids,
     * does k-means until converge, and prints results.
     *
     * @param argv      data file path and target number of clusters
     * @since           1.0
     */
    public static void main(String[] argv) throws Exception {
        
        if(argv.length != 2)
        {
            System.out.println("Proper Usage is: java ClusterDNA datafile #cluster");
            System.exit(0);
        }
        if (!(new File(argv[0]).isFile())) {
            System.out.println("Please check input file path");
            System.exit(0);
        }
        try{
            int checkInt=Integer.parseInt(argv[1]);
            if(checkInt<1){
                System.out.println("Please provide correct number of target clusters");
                System.exit(0);
            }
        }catch(NumberFormatException nfe){
            System.out.println("Please provide correct number of target clusters");
            System.exit(0);
        }
        
        
        /* the first argument is path of the data file containing the DNA strands */
        String inputFile=argv[0];
        /* the second argument is the number of the clusters */
        int k_param=Integer.parseInt(argv[1]);
        
        /* scan input file and store the DNA strands in a ArrayList*/
        ArrayList<ArrayList<String>> DNAs = new ArrayList<ArrayList<String>>();
        Scanner aScanner = new Scanner(new File(inputFile));
        String aLine = null;
        do {
            aLine = aScanner.nextLine();
            String[] parts = aLine.split(",");
            ArrayList<String> aDNA=new ArrayList<String>();
            for (int i=0;i<parts.length;i++){
                aDNA.add(parts[i]);
            }
            DNAs.add(aDNA);
        } while (aScanner.hasNext());
        /* get the length of each DNA strand */
        int lenDNA=DNAs.get(0).size();
        
        
        /* choose initial centroids, use either one of the two strategies below */
        
        /**
         * strategy 1
         * randomly select initial centroids from the original DNA strands,
         * and make sure the centroids are far away from each other
         */
        int maxSim=(int)(0.8f*lenDNA); // set it greater than maxSim in generatednadata.py
        ArrayList<ArrayList<String>> centroids = new ArrayList<ArrayList<String>>();
        int randInd=(int)(Math.random() * DNAs.size());
        centroids.add(DNAs.get(randInd));
        for (int i=0;i<k_param;i++){
            while(tooSimilar(DNAs.get(randInd),centroids,maxSim)){
                randInd=(int)(Math.random() * DNAs.size());
            }
            centroids.add(DNAs.get(randInd));
        }
        
        /**
         * strategy 2
         * use the first k DNA strands as initial centroids
         */
        /*
        ArrayList<ArrayList<String>> centroids = new ArrayList<ArrayList<String>>();
        for (int i=0;i<k_param;i++){
            centroids.add(DNAs.get(i));
        }
        */
        
        /* iterate until centroids are stable */
        boolean converged=false;
        while(!converged){
            
            /**
             * record number of different bases, 
             * in order to find max number of occurances in each location
             * and then determine the new centroid
             */
            int[][] tagCountA = new int[k_param][lenDNA];
            int[][] tagCountC = new int[k_param][lenDNA];
            int[][] tagCountG = new int[k_param][lenDNA];
            int[][] tagCountT = new int[k_param][lenDNA];
            
            /* record tag of each DNA strand */
            ArrayList<Integer> tags=new ArrayList<Integer>();
            
            /* assign tag for each DNA strand according to maximum similarity among centroids */
            for (int i=0;i<DNAs.size();i++){
                ArrayList<Integer> sims=new ArrayList<Integer>();
                for (int j=0;j<k_param;j++){
                    sims.add(simDNA(DNAs.get(i),centroids.get(j)));
                }
                int tag=sims.indexOf(Collections.max(sims));
                tags.add(tag);
                for(int j=0;j<lenDNA;j++){
                    if(DNAs.get(i).get(j).equals("A")){
                        tagCountA[tag][j] += 1;
                    }else if(DNAs.get(i).get(j).equals("C")){
                        tagCountC[tag][j] += 1;
                    }else if(DNAs.get(i).get(j).equals("G")){
                        tagCountG[tag][j] += 1;
                    }else if(DNAs.get(i).get(j).equals("T")){
                        tagCountT[tag][j] += 1;
                    }
                }
            }
            
            /* if the centroids are the same to the previous iteration, set converged */
            converged=true;
            for(int i=0;i<k_param;i++){
                ArrayList<String> updatedCentroid=new ArrayList<String>();
                for(int j=0;j<lenDNA;j++){
                    if(tagCountA[i][j]>=tagCountC[i][j] && tagCountA[i][j]>=tagCountG[i][j] && tagCountA[i][j]>=tagCountT[i][j]){
                        updatedCentroid.add("A");
                    }else if(tagCountC[i][j]>=tagCountA[i][j] && tagCountC[i][j]>=tagCountG[i][j] && tagCountC[i][j]>=tagCountT[i][j]){
                        updatedCentroid.add("C");
                    }else if(tagCountG[i][j]>=tagCountA[i][j] && tagCountG[i][j]>=tagCountC[i][j] && tagCountG[i][j]>=tagCountT[i][j]){
                        updatedCentroid.add("G");
                    }else if(tagCountT[i][j]>=tagCountA[i][j] && tagCountT[i][j]>=tagCountC[i][j] && tagCountT[i][j]>=tagCountG[i][j]){
                        updatedCentroid.add("T");
                    }
                    if(! updatedCentroid.get(j).equals(centroids.get(i).get(j))){
                        converged=false;
                    }
                }
                centroids.set(i,updatedCentroid);
            }
            
            /* if converged, print out the results */
            if(converged==true){
                for (int i=0;i<k_param;i++){
                    System.out.println("centroid "+i+":");
                    System.out.println(centroids.get(i));
                    System.out.println("points:");
                    for (int j=0;j<DNAs.size();j++){
                        if(tags.get(j)==i){
                            System.out.println(DNAs.get(j));
                        }
                    }
                }
            }
            
        }
        
    }
    
    
    /**
     * compute similarity between two DNA strands
     * 
     * @param p1    the DNA strand
     * @param p2    the other DNA strand
     * @return      the similarity between the two DNA strands
     * @since       1.0
     */
    public static Integer simDNA(ArrayList<String> p1, ArrayList<String> p2) {
        int dlen=p1.size();
        int sim=dlen;
        for(int ii=0;ii<dlen;ii++){
            if(!p1.get(ii).equals(p2.get(ii))){
                sim--;
            }
        }
        return sim;
    }
    
    
    /** 
     * decide if one DNA strand is too similar to any of the DNA strands in a list
     * 
     * @param one       the DNA strand
     * @param list      a list of DNA strand to be compared with
     * @param maxSim    the max number of identical bases allowed
     * @return          true if the DNA strand is too similar to one of the DNA strands in list
     * @since           1.0
     */
    public static boolean tooSimilar(ArrayList<String> one, ArrayList<ArrayList<String>> list, Integer maxSim){
        for (int iter=0;iter<list.size();iter++){
            if(simDNA(one,list.get(iter)) > maxSim){
                return true;
            }
        }
        return false;
    }
    
}

