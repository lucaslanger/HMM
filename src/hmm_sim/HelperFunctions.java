package hmm_sim;

import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class HelperFunctions {
	
	public static Random random = new Random();

	
	public static double[] generateProbabilityVector(int length, double p){	
		double[] result = new double[length]; 
		for(int i=0;i<length;i++){
			result[i] = nChooseI(length-1, i)*Math.pow(p,i)*Math.pow(1-p, length-1-i);
		} 
	
		return result;
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
	
	public static double[] listToArray(List<Double> arr){   
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
}
