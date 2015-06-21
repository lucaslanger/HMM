package hmm_sim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.Set;

import Jama.Matrix;

public class baumWelchExperimentation {
	
	private static Random random = new Random();
	
	//Development halted because EM not needed for project
	//TODO if needed: Verify correctness and add multiple observations
		
	public static void baumWelch(int numStates, int numIterations, int[] observations){
		double[][] prior =  { baumWelchExperimentation.randomNormVector(numStates) };
		double[][] transition = new double[numStates][numStates]; 
		double[][] observ = new double[][]{ baumWelchExperimentation.randomNormVector(2) , baumWelchExperimentation.randomNormVector(2) };
		
		for (int i = 0; i < numStates; i++) {
			transition[i] = baumWelchExperimentation.randomNormVector(numStates);
		}
			
		Matrix P,T,O,E;
		int duration = observations.length;
		double[][] forwardData = new double[ duration ][numStates];
		
		double[][][] transProbs = new double[duration-1][numStates][numStates];
		double[][] stateProbs = new double[duration][numStates];
		
		for (int c = 0; c < numIterations; c++) {
	
			P = new Matrix(prior);
			T = new Matrix(transition).transpose(); 
			O = new Matrix(observ);
			
			System.out.println("Run");
			P.print(numStates, numStates);
			T.print(numStates, numStates);
			O.print(numStates, numStates);
			
			int obs = 0;
			for (int i = 0; i < numStates; i++) {
				obs = observations[i];
				forwardData[0][i] = observ[obs][i] * prior[0][i];	
			}
			
			double tmp = 0;
			double tmp2 = 0;
			Matrix forw,tran;
			for (int i = 1; i < duration; i++) {
				for (int j = 0; j < numStates; j++) {
					tmp = 0;
					forw = new Matrix(forwardData[i-1],1).transpose();
					tran = new Matrix(transition[j], 1); 
					tmp = forw.times(tran).trace();
					
					tmp2 = 0;
					for (int p = 0; p < numStates; p++) {
						tmp += forwardData[i-1][p]*transition[j][p];	// from p to j
					}
					forwardData[i][j] = tmp*observ[obs][j];
				}
			}
			
			//System.out.println("Forward");
			//new Matrix(forwardData).print(5, 5);
			
			double[][] backwardData = new double[duration][numStates];
			for (int i = 0; i < numStates; i++) {
				backwardData[duration-1][i] = 1;
			}
			
			//Improvement: Turn the code below into matrix operations, numStates, observ states low so it doesn't matter
			int j;
			for (int i = 1; i < duration; i++) {
				for (int k = 0; k < numStates; k++) {
					tmp = 0;
					j = backwardData.length-1-i;
					for (int p = 0; p < numStates; p++) {
						if (observations[i] == 1){
							tmp += transition[p][k] * backwardData[j+1][p] * observ[0][p];
						}
						else{
							tmp += transition[p][k] * backwardData[j+1][p] * (1-observ[0][p]);
						}
					}
					if (tmp > 1) {
						//System.out.println(tmp);
					}
					
					backwardData[j][k] = tmp;
				}
			}
			
			//System.out.println("Backward");
			//new Matrix(backwardData).print(5, 5);
			
			double sum;
			for (int i = 0; i < duration; i++) {
				sum = 0;
				for (int k = 0; k < numStates; k++) {
					stateProbs[i][k] = forwardData[i][k]*backwardData[i][k];
					sum += stateProbs[i][k];
				}
				for (int k = 0; k < numStates; k++) {
					stateProbs[i][k] = stateProbs[i][k]/sum; 
				}
			}
			
			for (int i = 0; i < duration-1; i++) {
				sum = 0;
				for (int k = 0; k < numStates; k++) {
					for (int p = 0; p < numStates; p++) {
						if (observations[p] == 1){
							transProbs[i][k][p] = forwardData[i][k]*transition[p][k]*backwardData[i+1][p]*observ[0][p];
						}
						else{
							transProbs[i][k][p] = forwardData[i][k]*transition[p][k]*backwardData[i+1][p]*(1-observ[0][p]);
						}
						sum += transProbs[i][k][p];
					}
				}
				
				for (int k = 0; k < numStates; k++) {
					for (int p = 0; p < numStates; p++) {
						transProbs[i][k][p] = transProbs[i][k][p]/sum;
						if (transProbs[i][k][p] > 1) {
							System.out.println("Transprob");
							System.out.println(transProbs[i][k][p]);
						}
					}
				}
				
			}
			
			for (int i = 0; i < numStates; i++) {
				prior[0][i] = stateProbs[0][i];
			}
			
			double norm;
			for (int i = 0; i < numStates; i++) {
				for (int k = 0; k < numStates; k++) {
					norm = 0;
					for (int t = 0; t < duration-1; t++) {
						norm += stateProbs[t][i];
					}
					
					transition[i][k] = 0;
					for (int t= 0; t < duration-1; t++) {
						transition[i][k] += transProbs[t][i][k];
					}
					
					
					transition[i][k] = transition[i][k]/norm;
				}
			}
			
			double norm2, obscount;
			for (int i = 0; i < numStates; i++) {
				norm2 = 0;
				obscount = 0;
				for (int t = 0; t < duration; t++) {
					if( observations[t] == 1 ){
						obscount += stateProbs[t][i];
					}
					norm2 += stateProbs[t][i];
				}
				observ[0][i] = obscount/norm2;
			}
		}
		
		//Computation of P(Si, O1 .. Oi): Forward
		//Computation of P(Oi+1, ... On|Si) Backward
		//Now compute P(Si, O1...On) --> condition on Si, then you have a product
		//Transitions: P(Si,Si+1|O1...On), compute using info above and transitions guess
		//Special for Baum-Welch: since T's are fixed, compute sum over transitions from ij/ all possible transitions- average measure 
		//To compute O: taking sum t=1 to N #times Si=x*obs_frequency/ sum1toN obs_frequency
		//Idea above: Weighting observation likelyhood by likelyhood of Si=x there for t=1 to N
	}
	
	public static void testEMGaussian(){
		
		double c1 = random.nextDouble();
		double c2 = random.nextDouble();
		double c3 = random.nextDouble();
		
		System.out.println(c1);
		System.out.println(c2);
		System.out.println(c3);
		
		int N = 10;
		
		double[][] p = { {0}, {1}, {0} };
		double[][] t = { {0.3,0.3,0.4}, {0.3,0.3,0.4}, {0.3,0.3,0.4} }; 
		
		double[][] o = { 	baumWelchExperimentation.generateBinomialVector(N,c1), 
							baumWelchExperimentation.generateBinomialVector(N,c2), 
							baumWelchExperimentation.generateBinomialVector(N,c3) 
						};
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		
		HMM h = new HMM(T,O,P);
		
		int[] int_sequence = h.generateSequence(100);
		double[] sequence = new double[100];
		for (int i = 0; i < sequence.length; i++) {
			sequence[i] = (double) int_sequence[i];
		}
				
		Set<String> s = expectationMaximizationGaussian(sequence, 3 , 50);
		
		System.out.println( s );
	}
	
	public static Set<String> expectationMaximizationGaussian(double[] observationSequence, int clusters, int numIterations){
		double overallMean = getMean(observationSequence);
		double overallSD = getSd(observationSequence, overallMean);
		
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
				likelyhood = getLikelyhood(mean, sd, datapoint);
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
			l = doubleListToArray( emData.get(key));
			mean = getMean(l);
			sd = getSd(l, mean);
			String newkey = Double.toString(mean) + "," + Double.toString(sd);
			revisedEmData.put(newkey, new ArrayList<Double>());
		}
		return revisedEmData;
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
	
	public static double[] doubleListToArray(List<Double> arr){   
	    double[] result = new double[arr.size()];
	    int i = 0;
	    for(Double d : arr) {
	        result[i++] = d.doubleValue();
	    }
	    return result;
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
	
	public static double getLikelyhood(double modelMean, double modelSd, double obs){
		double num = Math.pow(Math.E, -1*Math.pow(modelMean - obs,2) / (2*Math.pow(modelSd,2)) );
		double denom = Math.sqrt(2*Math.pow(modelSd,2)*Math.PI);
		return num/denom;
	}
	
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
	
}
