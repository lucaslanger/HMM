package hmm_sim;

import java.util.HashMap;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class QueryEngine {
	
	private Matrix a0;
	private Matrix ainf;
	private Matrix[] Asigmas;
	private int maxExponent;
	private int maxPower;
	
	private boolean debug;
	//Fields below are only useful for debugging
	
	private Matrix pinv;
	private Matrix sinv;
	private HashMap<String, Matrix> truncatedSVD;
	private SingularValueDecomposition originalSVD;
	
	public QueryEngine(Matrix a0, Matrix ainf, Matrix[] Asigmas, int maxExponent, int base){
		this.a0 = a0;
		this.ainf = ainf;
		this.Asigmas = Asigmas;
		this.maxExponent = maxExponent;
		this.maxPower = (int) Math.pow(base, maxExponent);
		this.debug = false;
	}
	
	public QueryEngine(Matrix a0, Matrix ainf, Matrix[] Asigmas, int maxExponent, int base, Matrix pinv, Matrix sinv, HashMap<String, Matrix> truncatedSVD, SingularValueDecomposition originalSVD){
		this.a0 = a0;
		this.ainf = ainf;
		this.Asigmas = Asigmas;
		this.maxExponent = maxExponent;
		this.maxPower = (int) Math.pow(base, maxExponent);
		this.debug = true;
		
		this.pinv = pinv;
		this.sinv = sinv;
		this.truncatedSVD = truncatedSVD;
		this.originalSVD = originalSVD;
	}

	public Matrix getA0() {
		return a0;
	}

	public Matrix getAinf() {
		return ainf;
	}

	public Matrix[] getAsigmas() {
		return Asigmas;
	}

	public int getMaxExponent() {
		return maxExponent;
	}

	public int getMaxPower() {
		return maxPower;
	}
	
	public Matrix getH(){
		if (originalSVD == null) {
			System.out.println("Only non-debug QueryEngine instantiated");
			return null;
		}
		else{
			return originalSVD.getU().times(originalSVD.getS()).times(originalSVD.getV().transpose());
		}
	}
	
	/*public void debugHComparisons(){
		
		System.out.println("H differences");
		System.out.println(truncatedSVD.get(""));
		
		System.out.println("P's");
		machine.get("pinv").print(5, 5);
		tru.get("pinv").print(5, 5);
		machine.get("pinv").minus(tru.get("pinv")).print(5, 5);
		
		System.out.println("SVDs");
		machine.get("s_values").print(5, 5);
		tru.get("s_values").print(5, 5);	
		
		System.out.println("U's");
		machine.get("U").print(5, 5);
		tru.get("U").print(5, 5);
		
		System.out.println("VT's");

		machine.get("VT").print(5, 5);
		tru.get("VT").print(5, 5);
		
		System.out.println("Asigma=1 error");
		machine.get("1").minus(tru.get("1")).print(5,5);
	}*/

}
