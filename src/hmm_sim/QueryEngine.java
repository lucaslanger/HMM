package hmm_sim;

import java.util.HashMap;

import javax.management.Query;

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


	public static Matrix matrixPower(Matrix m, int exp){
		
		Matrix I = Matrix.identity(m.getArray().length, m.getArray().length);
		if (exp == 0){
			return I;
		}
		else{
			if (exp % 2 == 1) {
				return m.times( matrixPower(m, (exp-1)/2) );
			}
			else{
				return matrixPower(m, exp/2);
			}
		}
	}
	
	public Matrix matrixQuery(int power, int maxPower, int base, boolean forward){	
		int p = maxPower;
	
		int size = this.Asigmas[0].getArrayCopy().length;
		Matrix r = Matrix.identity(size,size);
		while(power != 0){
	
			while (p > power){
				p = p/base;
			}
			int exponent = (int)( Math.log(p)/Math.log(base));

			if (forward){
				r = r.times( this.Asigmas[exponent] );
			}
			else{
				r = this.Asigmas[exponent].times(r);
			}
			power -= p;
		}
		
		return r;
	}
	
	
	public Matrix alphaKQuery(int power, int maxPower, int base){	//Always forward for now
		
		int p = maxPower;
		
		Matrix r = this.a0;
		while(power != 0){
			
			while (p > power){
				p = p/base;
			}
		
			int exponent = (int)( Math.log(p)/Math.log(base));
		
			r = r.times( this.Asigmas[ exponent ] );
			power -= p;
		}
		
		return r;
		
	}

	public double probabilityQuery(int power, int maxPower, int base, boolean forward){
		int p = maxPower;
		Matrix r;
		if (forward){
			r = this.a0;
		}
		else{
			r = this.ainf;
		}
		while(power != 0){
			//r.print(5, 5);
			
			while (p > power){
				p = p/base;
			}
			int exponent = (int)( Math.log(p)/Math.log(base));
			
			if (forward){
				r = r.times( this.Asigmas[exponent] );
				
			}
			else{
				r = this.Asigmas[exponent].times(r);
			}
			power -= p;
		}
		if (forward){
			r = r.times(this.ainf);
		}
		else{
			r = this.a0.times(r);
		}
		try {
			QueryEngine.checkQueryResultValidity(r);
			return r.get(0,0);
		} catch (Exception e) {
			System.out.println("Problem with multiplication");
			e.printStackTrace();
			return -10000;
		}
		
		
	}
	
	private static void checkQueryResultValidity(Matrix r) throws Exception {

		if( !( (r.getArrayCopy().length == 1) && (r.getArrayCopy()[0].length == 1) ) ){
			r.print(5, 5);

			throw new Exception();
			
		}
		
	}

	public double debugProbabilityQuery(int power, int maxPower, int base, boolean forward){
		int p = maxPower;
		
		Matrix r;
		if (forward){
			r = this.a0;
		}
		else{
			r = ainf;
		}
		
		System.out.print("Exponent Ordering: ");
		while(power != 0){
			
			while (p > power){
				p = p/base;
			}
			int exponent = (int)( Math.log(p)/Math.log(base));
			
			System.out.print(p);
			System.out.print(", ");
			if (forward){	
				r = r.times( this.Asigmas[exponent] );	
			}
			else{
				r = this.Asigmas[exponent].times(r);
			}
			power -= p;
		}
		if (forward){
			r = r.times(this.ainf);
		}
		else{
			r = this.a0.times(r);
		}
		
		System.out.println("");
		
		return r.get(0,0);
		
	}
	
	public Matrix[] getAllKStateQueries(int maxK, int maxPow, int base){
		Matrix[] m = new Matrix[maxK];
		
		m[0] = this.a0;
		for (int i = 1; i < m.length; i++) {
			m[i] = this.alphaKQuery(i, maxPow, base );
			double normalizer = m[i].norm2();
			//double normalizer = m[i].times(ainf).get(0, 0);
			m[i] = m[i].times(1.0/normalizer);
			//System.out.println(m[i].norm2());
		}	
		return m;
	}
	
}
