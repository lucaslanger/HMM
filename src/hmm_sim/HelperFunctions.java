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
	
	public static double sumArray(double[] da){
		double s = 0;
		for(double d: da){
			s += d;
		}
		return s;
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
		
		int String = "Fix this slow procedure" + 111
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
		boolean debug = false;
		
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
	    
	    if (debug) {
	    	System.out.println("Before trunc");
	    	System.out.print("Size");
	    	System.out.println(S.getArrayCopy()[0].length);
	   	    //U.print(5, 5);
	   	    S.print(5, 5);
	   	    
	   	    //V.transpose().print(5, 5);
	   	    
	   	    System.out.println("After trunc");
	   	    System.out.print("Size");
	    	System.out.println(Strunc.getArrayCopy()[0].length);
	   	    //Utrunc.print(5, 5);
	   	    Strunc.print(5, 5);
	   	    //Vtrunc.print(5, 5);
		}
	 
	    
	    
	    return r;
	    
	}
	
	//TODO add paramter "maxdigit" to limit the query
	public static Matrix matrixQuery(HashMap<String, Matrix> d, int power, int maxPower, int base, boolean forward){	
		int p = maxPower;
		//System.out.println(d.get("max").get(0, 0));
		//System.out.println(d.keySet());
		int size = d.get( Integer.toString(p) ).getArrayCopy().length;
		Matrix r = Matrix.identity(size,size);
		while(power != 0){
	
			while (p > power){
				p = p/base;
			}
			if (forward){
				r = r.times( d.get(Integer.toString( p ) ) );
			}
			else{
				r = d.get(Integer.toString( p )).times(r);
			}
			power -= p;
		}
		
		return r;
	}
	
	public static double probabilityQuery(HashMap<String, Matrix> d, Matrix ao, Matrix ainf,  int power, int maxPower, int base, boolean forward){
		int p = maxPower;
		Matrix r;
		if (forward){
			r = ao;
		}
		else{
			r = ainf;
		}
		while(power != 0){
			
			while (p > power){
				p = p/base;
			}
			if (forward){
				//System.out.println(p);
				//System.out.println(d.keySet());
				r = r.times( d.get(Integer.toString( p ) ) );
				
			}
			else{
				r = d.get(Integer.toString( p )).times(r);
			}
			power -= p;
		}
		if (forward){
			r = r.times(ainf);
		}
		else{
			r = ao.times(r);
		}
		return r.get(0,0);
		
	}
	
	
	public static Matrix alphaKQuery(HashMap<String, Matrix> d, Matrix ao, int power, int maxPower, int base){	//Always forward for now
		
		int p = maxPower;
		
		Matrix r = ao;
		while(power != 0){
			
			while (p > power){
				p = p/base;
			}
		
			r = r.times( d.get(Integer.toString( p ) ) );
		
			power -= p;
		}
		
		return r;
		
	}
	
	public static double getMaxValue( double[] a){
		double max = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || max < a[i]){
				init = true;
				max = a[i];
			}
		}
		return max;
	}
	
	public static double getArgMax( double[] a){
		double max = 0;
		double argmax = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || max < a[i]){
				init = true;
				argmax = i;
				max = a[i];
			}
		}
		return argmax;
	}
	
	public static double getMinValue( double[] a){
		double min = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || min > a[i]){
				init = true;
				min = a[i];
			}
		}
		return min;
	}
	
	public static double getArgMin( double[] a){
		double min = 0;
		double argmin = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || min > a[i]){
				init = true;
				argmin = i;
				min = a[i];
			}
		}
		return argmin;
	}
	
	public static Matrix pseudoInvDiagonal(Matrix m){
		double[][] a = m.getArrayCopy();
		for (int i = 0; i < a.length; i++) {
			if (a[i][i] != 0){
				a[i][i] = 1/a[i][i];
			}
			else{
				a[i][i] = 0;
			}
		}
		return new Matrix(a);
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