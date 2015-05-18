//HMM model written by Lucas Langer to better my intuition about EM

package hmm_sim;

import Jama.*;

import java.util.*;

public class HMM {
	
	// T --> Transition, O --> Observation, P --> Prior 
	private Matrix T;
	private Matrix O;
	private Matrix P;
	
	public static Random random = new Random();
	
	public static void main(String[] args){
		
		
		double c1 = random.nextDouble();
		double c2 = random.nextDouble();
		double c3 = random.nextDouble();
		
		/*
		System.out.println(c1);
		System.out.println(c2);
		System.out.println(c3);
		*/
		
		double[][] t = { {0.3,0.3,0.4}, {0.3,0.3,0.4}, {0.3,0.3,0.4} }; 
		double[][] o = { generateBinomialVector(10,c1), 
						generateBinomialVector(10,c2), 
						generateBinomialVector(10,c3) };  
		
		double[][] p = { {0.3}, {0.4}, {0.4} };
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		
		HMM h = new HMM(T,O,P);
		
		double[] sequence = h.generateSequence(100);
		Set<String> s = h.expectationMaximizationGaussian(sequence, 3 , 50);
		
		testHankel();
	}

	public HMM(Matrix T, Matrix O, Matrix P){
		this.T = T;
		this.O = O;
		this.P = P;
	}
	
	public double[] generateSequence(int duration){
		double[] oSeq = new double[duration];
		
		int hiddenState = generateState( P.getArrayCopy()[0] );		
		
		for(int t=0;t<duration;t++){
			hiddenState = generateState( T.getArrayCopy()[hiddenState] );
			oSeq[t] = generateState( O.getArrayCopy()[hiddenState] );
		} 
		
		return oSeq;
	}

	public Set<String> expectationMaximizationGaussian(double[] observationSequence, int clusters, int numIterations){
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
	
	public HashMap<String, ArrayList<Double>> expect(HashMap<String, ArrayList<Double>> emData, double[] observationSequence){
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
	
	public HashMap<String, ArrayList<Double>> maximize(HashMap<String, ArrayList<Double>> emData){
		
		//Reason for copy is to avoid overwriting 
		HashMap<String, ArrayList<Double>> revisedEmData = new HashMap<String, ArrayList<Double>>();
		double[] l;
		double mean, sd;
		for (String key: emData.keySet()){
			l = listToArray( emData.get(key));
			mean = getMean(l);
			sd = getSd(l, mean);
			String newkey = Double.toString(mean) + "," + Double.toString(sd);
			revisedEmData.put(newkey, new ArrayList<Double>());
		}
		return revisedEmData;
	}
	
	public Matrix[] customEstimate(int numStates, Set<String> hiddenStateModels){
		//Use knowledge about data to infer on prior and then compute transition's for data manually
		return null;
	}
	
	public Matrix[] baumWelch(int numStates, int[] counts){
		//Seem's like a completely different algorithm because here we have no knowledge of model
		double[][] prior = new double[1][numStates];
		prior[0] = randomVector(numStates);
		
		double[][] transition = new double[numStates][numStates]; 
		double[][] observ = new double[1][numStates]; //Or continuous distribution
		observ[0] = randomVector(numStates);
		
		for (int i = 0; i < numStates; i++) {
			transition[i] = randomVector(numStates);
			
		}
		Matrix pi = new Matrix(prior);
		Matrix T = new Matrix(transition); 
		Matrix O = new Matrix(observ);
		
		int duration = counts[0];
		double[][] forwardData = new double[ duration ][numStates];
		for (int i = 0; i < forwardData.length; i++) {
			forwardData[0][i] = observ[0][i] * prior[0][i];	
		}
		
		double tmp = 0;
		for (int i = 1; i < duration; i++) {
			for (int j = 0; j < numStates; j++) {
				for (int p = 0; p < numStates; p++) {
					tmp += forwardData[i-1][p]*transition[p][j];
				}	
				forwardData[i][j] = tmp*observ[0][j];
				tmp = 0;
			}
		}
		
		double[][] backwardData = new double[duration][numStates];
		for (int i = 0; i < numStates; i++) {
			backwardData[duration-1][i] = 1;
		}
		
		//Improvement: Turn the code below into matrix operations, numStates, obs low it doesn't matter
		int j;
		for (int i = 1; i < duration; i++) {
			for (int k = 0; k < numStates; k++) {
				j = backwardData.length - i;
				for (int p = 0; p < numStates; p++) {
					tmp += transition[k][p] * backwardData[j+1][p] * observ[0][p];
				}
				backwardData[j][k] = tmp;
			}
		}
		
		double[][] temp1 = new double[duration][numStates];
		for (int i = 0; i < duration; i++) {
			double sum = 0;
			for (int k = 0; k < numStates; k++) {
				temp1[i][k] = forwardData[i][k]*backwardData[i][k];
				sum += temp1[i][k];
			}
			for (int k = 0; k < numStates; k++) {
				temp1[i][k] = temp1[i][k]/sum; 
			}
		}
		
		double[][][] temp2 = new double[duration-1][numStates][numStates];
		
		
		//Computation of P(Si, O1 .. Oi): Forward
		//Computation of P(Oi+1, ... On|Si) Backward
		//Now compute P(Si, O1...On) --> condition on Si, then you have a product
		//Transitions: P(Si,Si+1|O1...On), compute using info above and transitions guess
		//Special for Baum-Welch: since T's are fixed, compute sum over transitions from ij/ all possible transitions- average measure 
		//To compute O: taking sum t=1 to N #times Si=x*obs_frequency/ sum1toN obs_frequency
		//Idea above: Weighting observation likelyhood by likelyhood of Si=x there for t=1 to N
		
		return null;
	}

	
	public ArrayList<Matrix> singleObservationHankel(int[] counts, int base, int numHiddenStates){
	
		Matrix H = buildHankel(counts, 0);
		
		ArrayList<Matrix> H_Matrices  = new ArrayList<Matrix>();
		ArrayList<Matrix> A_Matrices  = new ArrayList<Matrix>();
		
		int maxDigit = (int) Math.floor( Math.log(counts.length )/Math.log(base) );
		int freq;
		for (int l = 1; l < maxDigit; l++) {
			freq = (int) Math.pow(base,l);
			H_Matrices.add( buildHankel(counts, freq) );
		}
		
		SingularValueDecomposition SVD = H.svd();
		Matrix pinv = SVD.getU().times(SVD.getS()).inverse();
		Matrix sinv = SVD.getV().inverse();
		
		for (int i = 0; i < H_Matrices.size(); i++) {
			A_Matrices.add( pinv.times(H_Matrices.get(i)).times( sinv ) );
		}
		
		double[][] h_L = new double[counts.length][1];
		for (int i = 0; i < counts.length; i++) {
			h_L[i][0] = (double) counts[i];
		}
		Matrix h_LS = new Matrix( h_L );
		Matrix h_PL = h_LS.transpose();
		
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_PL);
		
		A_Matrices.add(alpha_0);
		A_Matrices.add(alpha_inf);
	
		return A_Matrices;
	}

	public Matrix buildHankel(int[] counts, int startingindex, int endingindex){
		int size = endingindex - startingindex;
		double[][] hankel = new double[size][size];
		
		int i,j,k;
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				k = size - j;		
				hankel[k][j] = counts[i];
			}
		}
		
		return new Matrix(hankel);
	}
	
	public static void testHankel(){
		Matrix Hbar = new Matrix( new double[][]{ {0,0.2,0.14}, {0.2,0.22,0.15}, {0.14,0.45,0.31} }).transpose();
		
		Matrix Ha = new Matrix(new double[][]{ {0.2,0.22,0.15},{0.22,0.19,0.13},{0.45,0.45,0.32} }).transpose();
		Matrix Hb = new Matrix(new double[][]{ {0.14,0.45,0.31}, {0.15,0.29,0.13}, {0.31,0.85,0.58} } ).transpose();
		Matrix hls = new Matrix(new double[][]{ {0, 0.2, 0.14} } );
		Matrix hpl = new Matrix(new double[][]{ {0, 0.2, 0.14} } ).transpose();
		
		
		Hbar.print(5,5);
		/*
		Ha.print(5,5);
		Hb.print(5,5);
		*/
		
		SingularValueDecomposition svd = Hbar.svd();
		Matrix p = svd.getU().times( svd.getS() );
		Matrix s = svd.getV().transpose();
		
		Matrix pinv = p.inverse();
		Matrix sinv = s.inverse();
		
		Matrix Aa = pinv.times(Ha).times(sinv); 
		Matrix Ab = pinv.times(Hb).times(sinv); 
		
		Matrix alpha0 = hls.times(sinv);	//alpha0 row
		Matrix alphainf = pinv.times(hpl);	//alphainf column
		
		/*
		Aa.print(5,5);
		Ab.print(5,5);
		alpha0.print(5,5);
		alphainf.print(5,5);
		*/
		
		Matrix test1 = alpha0.times(Aa).times(alphainf);
		Matrix test2 = alpha0.times(Ab).times(alphainf);
		Matrix test3 = alpha0.times(Aa).times(Ab).times(alphainf);
		Matrix test4 = alpha0.times(Ab).times(Aa).times(alphainf);
		Matrix test5 = alpha0.times(Aa).times(Aa).times(alphainf);
		Matrix test6 = alpha0.times(Ab).times(Ab).times(alphainf);
		
		test1.print(5,5);
		test2.print(5,5);
		test3.print(5,5);
		test4.print(5,5);
		test5.print(5,5);
		test6.print(5,5);
	}
	
	
	
	// HELPER FUNCTIONS BELOW
	
	
	//Binomial Distribution
	public static double[] generateBinomialVector(int length, double p){	
		double[] result = new double[length]; 
		for(int i=0;i<length;i++){
			result[i] = nChooseI(length-1, i)*Math.pow(p,i)*Math.pow(1-p, length-1-i);
		} 
	
		return result;
	}
	
	//Uniform Random Probability Vector
	public static double[] randomVector(int size){	
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
