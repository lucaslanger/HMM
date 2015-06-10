package hmm_sim;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
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
	
	public static double[] incArray(int length){
		double[] output = new double[length];
		for (int i = 0; i < output.length; i++) {
			output[i] = i;
		}
		return output;
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

	
	public static HashMap<String, Matrix> truncateSVD(Matrix H, int nStates){
		
		SingularValueDecomposition svd = H.svd();
	    Matrix U = svd.getU();
	    Matrix S = svd.getS();
	    Matrix V = svd.getV();
	    
	    double[][] utemp = U.getArrayCopy();
	    double[][] utrunc = new double[utemp.length][nStates];
	    for (int i = 0; i < utrunc.length; i++) {
	    	for (int j = 0; j < nStates; j++) {
	    		utrunc[i][j] = utemp[i][j];
			}
			
		}
	    Matrix Utrunc = new Matrix(utrunc);
	    
	    double[][] stemp = S.getArrayCopy();
	    double[][] strunc = new double[nStates][nStates];
	    for (int i = 0; i < nStates; i++) {
	    	for (int j = 0; j < nStates; j++) {
	    		strunc[i][j] = stemp[i][j];
			}
		}
	    Matrix Strunc = new Matrix(strunc);
	    
	    double[][] vtemp = V.getArrayCopy();			//Double check to make sure this isnt going wrong
	    double[][] vtrunc = new double[utemp.length][nStates];
	    for (int i = 0; i < vtrunc.length; i++) {
	    	for (int j = 0; j < nStates; j++) {
	    		vtrunc[i][j] = vtemp[i][j];
			}
		}
	    Matrix Vtrunc = new Matrix(vtrunc).transpose();

	    HashMap<String, Matrix> r = new HashMap<String, Matrix>();
	    r.put("U", Utrunc);
	    r.put("VT", Vtrunc);
	    r.put("S", Strunc);
	    
	    /*
	    System.out.println("Before trunc");
	    U.print(5, 5);
	    S.print(5, 5);
	    V.transpose().print(5, 5);
	    
	    System.out.println("After trunc");

	    Utrunc.print(5, 5);
	    Strunc.print(5, 5);
	    Vtrunc.print(5, 5);
	    */
	    
	    return r;
	    
	}
	
	public static Matrix buildConcatHankel(double counts[], int basisSize){
		int l = counts.length-2*basisSize;
		double[][] concatHankel = new double[basisSize][basisSize*l];
		
		int h,s;
		for (int i = 0; i < basisSize; i++) {
			for (int j = 0; j < basisSize*l; j++) {
				h = j/basisSize;
				s = j%basisSize;
				/*System.out.println(i);
				System.out.println(j);
				System.out.println(basisSize);
				System.out.println(basisSize*l);
				System.out.println(h);
				System.out.println(s); */
				concatHankel[i][j] = counts[i + h + s ];
			}
		}
		return new Matrix(concatHankel);
		
	}
	
	//TODO add paramter "maxdigit" to limit the query
	public static Matrix matrixQuery(HashMap<String, Matrix> d, int power, int base, boolean forward){	
		int c = (int) d.get("max").norm1()-1;
		int p = (int) Math.pow(base, c);
		int size = d.get( Integer.toString(p) ).getArrayCopy().length;
		Matrix r = Matrix.identity(size,size);
		while(power != 0){
	
			while (p > power){
				p = p/base;
			}
			if (forward){
				r = r.times( d.get(Integer.toString( p ) ) );
				//System.out.println("Forward" );
			}
			else{
				r = d.get(Integer.toString( p )).times(r);
				//System.out.println("Backward" );
			}
			power -= p;
		}
		
		return r;
	}
	
	public static void outputData(String filename, String xaxisLabel, String yaxisLabel, double[][] xaxis, double[][] yaxis){
		try {
			PrintWriter writer = new PrintWriter(filename, "UTF-8");
			
			writer.println(xaxisLabel + "," + yaxisLabel);
			
			StringBuilder line = new StringBuilder();
			for (int j = 0; j < xaxis[0].length; j++) {	
				for (int i = 0; i < xaxis.length; i++) {
					line.append(xaxis[i][j] + ",");
					line.append(yaxis[i][j] + " ");
				}
				writer.println( line );
				line.setLength(0);
			} 
			
			writer.close();
			
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	


}