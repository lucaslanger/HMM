package hmm_sim;

import java.util.ArrayList;
import java.util.HashMap;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModel {
	
	private double[] probabilities;
	private SingularValueDecomposition svd;
	private int basisSize;

	
	public HankelSVDModel(double[] probabilities , int basisSize){
		this.probabilities = probabilities;
		this.basisSize = basisSize;
		this.svd = takeSVD();
	}
	
	public SingularValueDecomposition takeSVD(){
		try{
			Matrix H = buildH(0, basisSize);
			SingularValueDecomposition svd = H.svd();
			return svd;
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}
	
	public HashMap<String, Matrix> buildHankelBasedModel(int basisSize, int base, int modelSize){
		
		HashMap<String, Matrix> truncatedSVD = this.truncateSVD(modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
						
		ArrayList<Matrix> H_Matrices = new ArrayList<Matrix>();
		HashMap<String, Matrix> modelInformation = new HashMap<String, Matrix>();
		
		int maxDigit = (int) Math.floor((Math.log( (this.probabilities.length/2) - basisSize)/Math.log(base))) ; 
		int freq;
		Matrix h;
		for (int l = 0; l <= maxDigit; l++) {
			freq = (int) Math.pow(base,l);
			try {
				h = this.buildH(freq, freq+basisSize);
				H_Matrices.add(h);
			} catch (Exception e) {
				System.out.println("Problem Building Model when creating Hankel");
				e.printStackTrace();
				return null;
			}
		}
		
		Matrix AsigmaX;
		for (int i = 0; i < H_Matrices.size(); i++) {
			AsigmaX = pinv.times(H_Matrices.get(i)).times( sinv );
			
			int X = (int) Math.pow(2, i);
			modelInformation.put( Integer.toString(X), AsigmaX);
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < basisSize; i++) {
			h_L[i][0] = this.probabilities[i];
		}
		
		Matrix h_LS = new Matrix( h_L ).transpose();
		Matrix h_PL = h_LS.transpose();
				
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_PL);
				
		Matrix maxExponent = new Matrix(new double[][]{{maxDigit}});
		
		modelInformation.put("maxExponent", maxExponent);
		modelInformation.put("a0", alpha_0);
		modelInformation.put("ainf", alpha_inf);
		
		/*	//Only for debugging
		 * 		
		Matrix hTruncated = truncatedSVD.get("U").times(truncatedSVD.get("S")).times(truncatedSVD.get("VT"));
		modelInformation.put("Htruncated", hTruncated);
		modelInformation.put("pinv", pinv);
		modelInformation.put("sinv", sinv);
		modelInformation.put("S", truncatedSVD.get("S"));
		modelInformation.put("U", truncatedSVD.get("U"));
		modelInformation.put("VT", truncatedSVD.get("VT"));
		*/

		return modelInformation;
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
	
	public HashMap<String, Matrix> truncateSVD(int nStates){
		boolean debug = false;
		
	    Matrix U = this.svd.getU();
	    Matrix S = this.svd.getS();
	    Matrix V = this.svd.getV();
	    
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
		

}
