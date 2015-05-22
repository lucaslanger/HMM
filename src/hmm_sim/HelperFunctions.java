package hmm_sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Set;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HelperFunctions {
	
	public static Random random = new Random();	
	
	public static Set<String> expectationMaximizationGaussian(double[] observationSequence, int clusters, int numIterations){
		double overallMean = HelperFunctions.getMean(observationSequence);
		double overallSD = HelperFunctions.getSd(observationSequence, overallMean);
		
		HashMap<String, ArrayList<Double>> emData = new HashMap<String, ArrayList<Double>>();
		
		double meanGuess, sdGuess;
		for(int count = 0;count < clusters;count++){
		    meanGuess = overallMean - overallSD + random.nextDouble()*4*overallSD;	//4 arbitary here
			sdGuess = overallSD;
			
			String key = Double.toString(meanGuess) + "," + Double.toString(sdGuess);
			emData.put(key, new ArrayList<Double>() ); 
		}
		
		for(int iteration=0;iteration<numIterations;iteration++){
			
			expect(emData, observationSequence);
			//Investigate: WHY is emData not mutable under maximize?
			emData = maximize(emData);
		}
		
		return emData.keySet();
	}
	
	public static HashMap<String, ArrayList<Double>> expect(HashMap<String, ArrayList<Double>> emData, double[] observationSequence){
		for (String k: emData.keySet()){
			emData.put(k, new ArrayList<Double>());
		}
		
		String maxKey = null;
		double mean, sd, likelyhood, maxProbability = 0;
		ArrayList<Double> currentData;
		for(double datapoint: observationSequence){
			for(String key: emData.keySet()){
				mean = Double.parseDouble(key.split(",")[0]);
				sd = Double.parseDouble(key.split(",")[1]);
				likelyhood = HelperFunctions.getLikelyhood(mean, sd, datapoint);
				if (maxKey == null || likelyhood > maxProbability ){
					maxKey = key;
					maxProbability = likelyhood;
				}
			}
			currentData = emData.get(maxKey);
			currentData.add(datapoint);
			
			emData.put(maxKey, currentData);
			
			maxProbability = 0;
			maxKey = null;
		}
		
		return emData;

	}
	
	public static HashMap<String, ArrayList<Double>> maximize(HashMap<String, ArrayList<Double>> emData){
		
		//Reason for copy is to avoid overwriting 
		HashMap<String, ArrayList<Double>> revisedEmData = new HashMap<String, ArrayList<Double>>();
		double[] l;
		double mean, sd;
		for (String key: emData.keySet()){
			l = HelperFunctions.doubleListToArray( emData.get(key));
			mean = HelperFunctions.getMean(l);
			sd = HelperFunctions.getSd(l, mean);
			String newkey = Double.toString(mean) + "," + Double.toString(sd);
			revisedEmData.put(newkey, new ArrayList<Double>());
		}
		return revisedEmData;
	}
	
	//From Binomial Distribution
	public static double[] generateBinomialVector(int length, double p){	
		double[] result = new double[length]; 
		for(int i=0;i<length;i++){
			result[i] = nChooseI(length-1, i)*Math.pow(p,i)*Math.pow(1-p, length-1-i);
		} 
	
		return result;
	}
	
	//Uniform Random Probability Vector
	public static double[] randomNormVector(int size){	
		double[] vector = new double[size];
		double sum = 0;
		for (int i = 0; i < vector.length; i++) {
			vector[i] = random.nextDouble();
			sum += vector[i];
		}
		for (int j = 0; j < vector.length; j++) {	//Normalize
			vector[j] = vector[j]/sum;
		}
		return vector;
		
	}
	
	public static double[] randomProbVector(int size){
		double[] r = new double[size];
		for (int i = 0; i < size; i++) {
			r[i] = random.nextDouble();
		}
		return r;
	}
	
	//[0.4, 0.5, 0.1] --> returns index in {0,1,2}
	public static int generateState(double[] stateProbabilities){ 						
		int l = stateProbabilities.length; 
		double[] cumulativeSum = new double[l];
		
		cumulativeSum[0] = stateProbabilities[0];
		for (int i = 1; i<l;i++){
			cumulativeSum[i] = cumulativeSum[i-1] + stateProbabilities[i];
		}

		double r = random.nextDouble();
		
		int index = Arrays.binarySearch(cumulativeSum, r);

		if (index >= 0){
			return index;
		}
		else{
			return -1*(index + 1);
		}
		
	}	
	
	public static double[] doubleListToArray(List<Double> arr){   
	    double[] result = new double[arr.size()];
	    int i = 0;
	    for(Double d : arr) {
	        result[i++] = d.doubleValue();
	    }
	    return result;
	}
	
	public static double sumArray(double[] da){
		double s = 0;
		for(double d: da){
			s += d;
		}
		return s;
	}
	
	
	//FIX ELSE STATEMENT to generalize
	public static double getMean(double[] da){
		double m = 0;
		for (double d: da){
			m += d;
		}
		if (da.length != 0){
			return m/da.length;
		}
		else{
			return random.nextDouble()*2;
		}
	}
	
	//FIX ELSE STATEMENT to generalize
	public static double getSd(double[] da, double mean){
		double sd = 0;
		for (double d: da){
			sd += Math.pow(d-mean, 2);
		}
		if (sd!=0){
			return Math.sqrt(sd/da.length);
		}
		else{	//Quickfix to avoid 0 sd
			return 0.01; 
		}
	}
	
	// Gaussian pdf
	public static double getLikelyhood(double modelMean, double modelSd, double obs){
		double num = Math.pow(Math.E, -1*Math.pow(modelMean - obs,2) / (2*Math.pow(modelSd,2)) );
		double denom = Math.sqrt(2*Math.pow(modelSd,2)*Math.PI);
		return num/denom;
	}
	
	public static int nChooseI(int n, int i){
		return factorial(n)/ ( factorial(i)*factorial(n-i) );
	}
	public static int factorial(int n){
		if (n == 0){
			return 1;
		}
		else{
			return n*factorial(n-1);
		}
	}
	
	public static Matrix matrixPower(Matrix m, int exp){
		
		
		Matrix I = Matrix.identity(m.getArray().length, m.getArray().length);
		Matrix T = I;
		if (exp == 0){
			return I;
		}
		else{
			for (int i = 0; i < exp; i++) {
				T = T.times(m);
			}
			return T;
		}
		/*
		else{
			SingularValueDecomposition s = m.svd();
		} */
	}
}
