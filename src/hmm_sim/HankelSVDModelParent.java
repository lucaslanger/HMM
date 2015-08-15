package hmm_sim;

import java.io.IOException;
import java.io.Serializable;
import java.util.HashMap;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public abstract class HankelSVDModelParent implements Serializable{
	
	private int basisSize;
	private SingularValueDecomposition svdOfH;
	
	
	public abstract QueryEngine buildHankelBasedModel(int base, int modelSize);
		
	
	public static Matrix pseudoInvDiagonalKL(Matrix m){
		double[][] a = m.getArrayCopy();
		for (int i = 0; i < a.length; i++) {
			if (a[i][i] > m.get(0,0)*0.001){
				a[i][i] = 1/a[i][i];
			}
			else{
				a[i][i] = 0;
			}
		}
		
		Matrix r = new Matrix(a);
		/*System.out.println("Old Rank");
		System.out.println(m.rank());
		System.out.println("New Rank");
		System.out.println(r.rank());
		System.out.println();
		*/
		return r;
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
		
		Matrix r = new Matrix(a);
		
		if ( sumOfDiagonal(m.times(r)) > m.rank()*1.01 ){
			System.out.println("Weird numerical behavior in computing the pseudo diagonal inverse!");
			System.out.println(sumOfDiagonal(m.times(r)));
			System.out.println(r);
			System.out.println();
		}
		/*
		else{
			System.out.println("Working just fine");
			System.out.println(sumOfDiagonal(m.times(r)));
			System.out.println(m.rank());
			System.out.println();
		}
		*/
		return r;
	}
	
	public static double sumOfDiagonal(Matrix m){
		double sum = 0;
		for (int i = 0; i < m.getArrayCopy().length; i++) {
			sum += m.get(i, i);
		}
		return sum;
	}


	public static double[] getDiagonalArray(Matrix m){
		double[] r = new double[m.getArrayCopy()[0].length];
		for (int i = 0; i < r.length; i++) {
			r[i] = m.get(i, i);
		}
		return r;
	}
	
	public int getRank(){
		return this.svdOfH.getS().rank();
	}
	
	public Matrix getHankel(){
		return this.svdOfH.getU().times(this.svdOfH.getS()).times(this.svdOfH.getV().transpose());
	}

	public static HashMap<String, Matrix> takeSVD(Matrix H){
		//printMatrixDimensions(H, "Hankel");
		
		HashMap<String, Matrix> s = new HashMap<String, Matrix>();
		if (H.getArrayCopy().length <= H.getArrayCopy()[0].length){
			//System.out.println("Doing tranpose trick! \n");
			SingularValueDecomposition svd = new SingularValueDecomposition(H.transpose());
			
			//System.out.println("Careful with assignments below, make sure they are correct!");
			//System.out.println();
			s.put("S", svd.getS().transpose());
			s.put("U", svd.getV());
			s.put("VT", svd.getU().transpose());
			
			/*printMatrixDimensions(s.get("U"), "U");
			printMatrixDimensions(s.get("S"), "S");
			printMatrixDimensions(s.get("VT"), "VT");
			*/
		}
		else{
			SingularValueDecomposition svd = new SingularValueDecomposition(H);
			s.put("S", svd.getS());
			s.put("U", svd.getU());
			s.put("VT", svd.getV().transpose());
			
		}
		return s;
	}
	
	private static void printMatrixDimensions(Matrix h, String id){
		//System.out.println("Number of columns: " + h.getArrayCopy()[0].length);
		//System.out.println("Number of rows: " + h.getArrayCopy().length);
		System.out.println( "Dimensions of: " + id + " " + h.getArrayCopy().length + "x" + h.getArrayCopy()[0].length);
		System.out.println();
	}
	
	public static HashMap<String, Matrix> truncateSVD(Matrix H, int nStates){
		
		HashMap<String, Matrix> svdLocal = takeSVD(H);
	    Matrix U = svdLocal.get("U");
	    Matrix S = svdLocal.get("S");
	    Matrix V = svdLocal.get("VT").transpose();
	    
	    /*System.out.println("TESTING TRUNCATION");
	    printMatrixDimensions(H, "H");
	    printMatrixDimensions(U, "U");
	    printMatrixDimensions(S, "S");
	    printMatrixDimensions(V.transpose(), "VT");
	  	*/
	  
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
	    double[][] vtrunc = new double[vtemp.length][nStates];
	    for (int i = 0; i < vtrunc.length; i++) {
	    	for (int j = 0; j < nStates; j++) {
	    		vtrunc[i][j] = vtemp[i][j];
			}
		}
	    Matrix VTtrunc = new Matrix(vtrunc).transpose();

	    /*System.out.println("DIVSION");
	    SingularValueDecomposition svdToReturn = new SingularValueDecomposition(Utrunc.times(Strunc).times(VTtrunc));
	    System.out.println("TESTING HANKELSVDMODEPARENT");
	    System.out.println(Utrunc.rank());
	    System.out.println(svdToReturn.getU().rank());
	    System.out.println(Utrunc.minus(svdToReturn.getU()).norm2());
	    System.out.println(Strunc.minus(svdToReturn.getS()).norm2());
	    System.out.println(VTtrunc.minus(svdToReturn.getV().transpose()).norm2());

	    */
	    HashMap<String, Matrix> r = new HashMap<String, Matrix>();
	    r.put("U", Utrunc);
	    r.put("S", Strunc);
	    r.put("VT", VTtrunc);
	    
	    /*printMatrixDimensions(Utrunc, "Utrunc");
	    printMatrixDimensions(Strunc, "Strunc");
	    printMatrixDimensions(VTtrunc, "VTtrunc");
	    */
	    
	    boolean debug = false;
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
}
