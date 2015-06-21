package hmm_sim;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.HashMap;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModel {
	
	private double[] probabilities;
	
	private Matrix Prefixes;
	private Matrix Suffixes;

	
	public HankelSVDModel(double[] probabilities){
		this.probabilities = probabilities; 
	}
	
	public HashMap<String, Matrix> buildHankelBasedModel(int basisSize, int base, int numHiddenStates){
		
		Matrix H = buildH(0, basisSize);
		
		HashMap<String, Matrix> SVD = HelperFunctions.truncateSVD(H, numHiddenStates);
		Matrix Htrunc = SVD.get("U").times(SVD.get("S")).times(SVD.get("VT"));	
	
		Matrix di = HelperFunctions.pseudoInvDiagonal(SVD.get("S"));
		Matrix pinv = di.times(SVD.get("U").transpose());
		Matrix sinv = (SVD.get("VT")).transpose();
						
		ArrayList<Matrix> H_Matrices  = new ArrayList<Matrix>();
		HashMap<String, Matrix> returnData  = new HashMap<String, Matrix>();
		
		int maxDigit = (int) Math.floor((Math.log( (counts.length/2) - basisSize)/Math.log(base))) ; //Too low fix later to allow higher powers
		int freq;
		Matrix h;
		for (int l = 0; l <= maxDigit; l++) {
			freq = (int) Math.pow(base,l);
			h = buildHankel(counts, freq, freq+basisSize);
			H_Matrices.add(h);
		}
		
		Matrix m;
		for (int i = 0; i < H_Matrices.size(); i++) {
			m = pinv.times(H_Matrices.get(i)).times( sinv );
			returnData.put(Integer.toString( (int) Math.pow(2, i) ), m );
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < basisSize; i++) {
			h_L[i][0] = (double) counts[i];
		}
		Matrix h_LS = new Matrix( h_L ).transpose();
		Matrix h_PL = h_LS.transpose();
				
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_PL);
				
		Matrix maX = new Matrix(new double[][]{{maxDigit}});
		
		returnData.put("H", Htrunc);
		returnData.put("max", maX);
		returnData.put("pinv", pinv);
		returnData.put("sinv", sinv);
		returnData.put("S", SVD.get("S"));
		returnData.put("U", SVD.get("U"));
		returnData.put("VT", SVD.get("VT"));
		returnData.put("a0", alpha_0);
		returnData.put("ainf", alpha_inf);
		
		Matrix dA = new Matrix(new double[][]{{dataAmount}});
		returnData.put("dataAmount", dA);

		return returnData;
	}
	
	public Matrix buildH(int startingIndex, int endingIndex) throws Exception{		
		int hSize = (endingIndex - startingIndex);
		if ( (hSize + startingIndex)*2 > this.probabilities.length){
			throw new Exception("You asked for too large a Hankel Matrix. Increase the max durations recorded in your Environment.");
		} 
		else{
			double[][] hankel = new double[hSize][hSize];
			for (int i = 0; i < hSize; i++) {
				for (int j = 0; j < hSize; j++) {
					hankel[i][j] = this.probabilities[i+j+startingIndex];
				}
			}
			Matrix H = new Matrix(hankel);
			return H;
		}
		
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
	
	public void getTruncatedModel(int nStates){
		
	}
	
	public static void testHankel(){
		Matrix Hbar = new Matrix( new double[][]{ {0,0.2,0.14}, {0.2,0.22,0.15}, {0.14,0.45,0.31} }).transpose();
		
		Matrix Ha = new Matrix(new double[][]{ {0.2,0.22,0.15},{0.22,0.19,0.13},{0.45,0.45,0.32} }).transpose();
		Matrix Hb = new Matrix(new double[][]{ {0.14,0.45,0.31}, {0.15,0.29,0.13}, {0.31,0.85,0.58} } ).transpose();
		Matrix hls = new Matrix(new double[][]{ {0, 0.2, 0.14} } );
		Matrix hpl = new Matrix(new double[][]{ {0, 0.2, 0.14} } ).transpose();
			
		Hbar.print(5,5);
		
		SingularValueDecomposition svd = Hbar.svd();
		Matrix p = svd.getU().times( svd.getS() );
		Matrix s = svd.getV().transpose();
		
		Matrix pinv = p.inverse();
		Matrix sinv = s.inverse();
		
		Matrix Aa = pinv.times(Ha).times(sinv); 
		Matrix Ab = pinv.times(Hb).times(sinv); 
		
		Matrix alpha0 = hls.times(sinv);	//alpha0 row
		Matrix alphainf = pinv.times(hpl);	//alphainf column
		
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
		

}
