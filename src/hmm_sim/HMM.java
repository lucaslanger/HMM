//HMM model written by Lucas Langer to better one's intuition about EM

package hmm_sim;

import Jama.*;

import java.util.*;

public class HMM {
	
	// T --> Transition, O --> Observation, P --> Prior 
	private Matrix T;
	private Matrix O;
	private Matrix P;
	
	private Random random;
	
	public HMM(Matrix T, Matrix P, Matrix O){
		this.T = T;
		this.O = O;
		this.P = P;
	}
	
	public double[] generateSequence(int duration){
		double[] oSeq = new double[duration];
		
		int hiddenState = generateState(this.P.getArray()[0]);		
		int t = 0;
		
		while (t < duration){
			hiddenState = generateState( T.getArray()[hiddenState] );
			oSeq[t] = generateState( O.getArray()[hiddenState] );
		} 
		
		return oSeq;
	}

	//[0.4, 0.5, 0.1] --> returns index in {0,1,2}
	public int generateState(double[] stateProbabilities){ 						
		int l = stateProbabilities.length; 
		double[] cumulativeSum = new double[l-1];
		
		cumulativeSum[0] = stateProbabilities[0];
		for (int i = 1; i<l-1;i++){
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
	
	
	//Move this to another class?
	public Matrix expectationMaximizationGaussian(double[] observationSequence, int clusters, int numIterations){
		double overallMean = getMean(observationSequence);
		double overallSD = getSd(observationSequence, overallMean);
		
		HashMap<String, ArrayList<Double>> emData = new HashMap<String, ArrayList<Double>>();
		
		for(int count = 0;count < clusters;count++){
			double meanGuess = overallMean - overallSD + random.nextDouble()*2*overallSD;
			double sdGuess = overallSD/clusters;
			
			String key = Double.toString(meanGuess) + "," + Double.toString(sdGuess);
			emData.put(key, new ArrayList<Double>() );  	
		}
		
		int iteration = 0;
		while (iteration < numIterations){
			
			emData = maximize(emData, observationSequence);
			expect(emData);
			
		}
		
		
	}
	
	public HashMap<String, ArrayList<Double>> expect(HashMap<String, ArrayList<Double>> emData, double[] observationSequence){
		
	}
	
	public HashMap<String, ArrayList<Double>> maximize(HashMap<String, ArrayList<Double>> emData, double[] observationSequence){
		
		HashMap<String, ArrayList<Double>> revisedEmData = new HashMap<String, ArrayList<Double>>();
		for (String key: emData.keySet()){
			double[] l = listToArray( emData.get(key));
			double mean = getMean(l);
			double sd = getSd(l, mean);
			String newkey = Double.toString(mean) + "," + Double.toString(sd);
			revisedEmData.put(newkey, new ArrayList<Double>());
		}
		
		return revisedEmData;
	}
	
	public static double[] listToArray(List<Double> arr){   
	    double[] result = new double[arr.size()];
	    int i = 0;
	    for(Double d : arr) {
	        result[i++] = d.doubleValue();
	    }
	    return result;
	}
	
	
	public double getMean(double[] da){
		double m = 0;
		for (double d: da){
			m += d;
		}
		return m/da.length;
	}
	
	public double getSd(double[] da, double mean){
		double m = 0;
		for (double d: da){
			m += Math.pow(d-mean, 2);
		}
		return Math.sqrt(m/da.length);
	}
	
	public double getLikelyhood(double modelMean, double modelSd, double obs){
		double num = Math.pow(Math.E, -1*Math.pow(modelMean - obs,2) / (2*Math.pow(modelSd,2)) );
		double denom = Math.sqrt(2*Math.pow(modelSd,2)*Math.PI);
		return num/denom;
	}

}
