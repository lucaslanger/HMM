package hmm_sim;

import java.io.IOException;
import java.util.HashMap;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public abstract class HankelSVDModelParent {
	
	private int basisSize;
	private SingularValueDecomposition svdOfH;
	
	
	public abstract QueryEngine buildHankelBasedModel(int base, int modelSize);
	
	public abstract SingularValueDecomposition truncateSVD(int nStates);
	
	
	public abstract synchronized void writeObject(java.io.ObjectOutputStream stream);
	
	public abstract void readObject(java.io.ObjectInputStream in);
	
	public abstract void test();
	
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

}
