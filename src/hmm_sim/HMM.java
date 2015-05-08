//HMM model written by Lucas Langer to better my intuition about EM

package hmm_sim;

import Jama.*;

import java.util.*;

public class HMM {
	
	// T --> Transition, O --> Observation, P --> Prior 
	private Matrix T;
	private Matrix O;
	private Matrix P;
	
	public Random random;
	
	public static void main(String[] args){
		
		
		double[][] t = { {0.6,0.3,0.1}, {0.3,0.3,0.4}, {0.1,0.7,0.2} }; 
		double[][] o = { {0.2,0.3,0.5}, {0.3,0.3,0.4}, {0.1,0.7,0.2} };  
		
		double[][] p = { {0.3}, {0.4}, {0.4} };
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		double i = O.getArrayCopy()[1][1];
		
		Matrix P = new Matrix( p );
		HMM h = new HMM(T,O,P);
		
		double[] sequence = h.generateSequence(100);
		System.out.println( Arrays.toString(sequence) );
		Set<String> s = h.expectationMaximizationGaussian(sequence, 3 , 50);
		System.out.println( s.toString() );
		
		
		//testHankel();
	}

	public HMM(Matrix T, Matrix O, Matrix P){
		this.T = T;
		this.O = O;
		this.P = P;
		
		random = new Random();
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

	//[0.4, 0.5, 0.1] --> returns index in {0,1,2}
	public int generateState(double[] stateProbabilities){ 						
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
	
	//Move this to another class?
	public Set<String> expectationMaximizationGaussian(double[] observationSequence, int clusters, int numIterations){
		double overallMean = getMean(observationSequence);
		double overallSD = getSd(observationSequence, overallMean);
		
		HashMap<String, ArrayList<Double>> emData = new HashMap<String, ArrayList<Double>>();
		
		double meanGuess, sdGuess;
		for(int count = 0;count < clusters;count++){
		    meanGuess = overallMean - overallSD + random.nextDouble()*2*overallSD;
			sdGuess = overallSD/clusters;
			
			String key = Double.toString(meanGuess) + "," + Double.toString(sdGuess);
			emData.put(key, new ArrayList<Double>() ); 
		}
		
		for(int iteration=0;iteration<numIterations;iteration++){
			expect(emData, observationSequence);
			emData = maximize(emData);
			System.out.println( emData.keySet() );
		}
		
		return emData.keySet();
	}
	
	public void expect(HashMap<String, ArrayList<Double>> emData, double[] observationSequence){
		String maxKey = null;
		double maxProbability = 0;
		
		double mean, sd, likelyhood;
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
		emData = revisedEmData;
		return emData;
	}
	
	public double getLikelyhood(double modelMean, double modelSd, double obs){
		double num = Math.pow(Math.E, -1*Math.pow(modelMean - obs,2) / (2*Math.pow(modelSd,2)) );
		double denom = Math.sqrt(2*Math.pow(modelSd,2)*Math.PI);
		return num/denom;
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
	
	
	public HashMap<String, Matrix> hankelSVD(HashMap<String, Double> fvals){
		//CREATE HANKEL CODE HERE: 
		//SingularValueDecomposition svd = hankelMatrix.svd();
		//Matrix p = svd.getU().times( svd.getS() );
		//Matrix s = svd.getV();
		return null;
	}
	
	public static void testHankel(){
		Matrix Hbar = new Matrix( new double[][]{ {0,0.2,0.14}, {0.2,0.22,0.15}, {0.14,0.45,0.31} }).transpose();
		
		Matrix Ha = new Matrix(new double[][]{ {0.2,0.22,0.15},{0.22,0.19,0.13},{0.45,0.45,0.32} }).transpose();
		Matrix Hb = new Matrix(new double[][]{ {0.14,0.45,0.31}, {0.15,0.29,0.13}, {0.31,0.85,0.58} } ).transpose();
		Matrix hls = new Matrix(new double[][]{ {0, 0.2, 0.14} } );
		Matrix hpl = new Matrix(new double[][]{ {0, 0.2, 0.14} } ).transpose();
		
		//hls.print(5,5);
		
		SingularValueDecomposition svd = Hbar.svd();
		Matrix p = svd.getU().times( svd.getS() );
		Matrix s = svd.getV();
		Matrix pinv = p.inverse();
		Matrix sinv = s.inverse();
		
		Matrix Aa = pinv.times(Ha).times(sinv); 
		Matrix Ab = pinv.times(Hb).times(sinv); 
		
		Matrix alpha0 = hls.times(sinv);	//alpha0 row
		Matrix alphainf = pinv.times(hpl);	//alphainf column
		
		Aa.print(5,5);
		Ab.print(5,5);
		alpha0.print(5,5);
		alphainf.print(5,5);
		
		Matrix r = alpha0.times(Ab).times(alphainf);
		r.print(5,5);
	}
	
	
}
