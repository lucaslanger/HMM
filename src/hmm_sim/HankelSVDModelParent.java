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

	
	public static HashMap<String, Matrix> truncateSVD(Matrix H, int nStates){
		
		SingularValueDecomposition svdLocal = H.svd();
	    Matrix U = svdLocal.getU();
	    Matrix S = svdLocal.getS();
	    Matrix V = svdLocal.getV();
	    
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
	    r.put("VT", VTtrunc);
	    r.put("S", Strunc);
	    
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
