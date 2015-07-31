package hmm_sim;


import java.util.HashMap;

import javax.management.Query;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class QueryEngineMultipleObservations {
		
	private Matrix a0;
	private Matrix ainf;
	private HashMap<String, Matrix> Asigmas;
	private int maxExponent;
	private int maxPower;
	
	private boolean debug;
	//Fields below are only useful for debugging
	
	private Matrix pinv;
	private Matrix sinv;
	private HashMap<String, Matrix> truncatedSVD;
	private SingularValueDecomposition originalSVD;
	
	public QueryEngineMultipleObservations(Matrix a0, Matrix ainf, HashMap<String, Matrix> Asigmas, int base){
		this.a0 = a0;
		this.ainf = ainf;
		this.Asigmas = Asigmas;

		this.debug = false;
	}
	
	public QueryEngineMultipleObservations(Matrix a0, Matrix ainf, HashMap<String, Matrix> Asigmas, int base, Matrix pinv, Matrix sinv, HashMap<String, Matrix> truncatedSVD, SingularValueDecomposition originalSVD){
		this.a0 = a0;
		this.ainf = ainf;
		this.Asigmas = Asigmas;
		
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

	public HashMap<String, Matrix> getAsigmas() {
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

	public double probabilityQuery(SequenceOfSymbols sequence){
		Matrix r = this.a0;
		while(sequence.getSequence() != ""){
			
			SequenceOfSymbols nextstreak = sequence.getFirstStreak();
			int power = nextstreak.getStreakFromString();
			String symbol = nextstreak.getSymbolFromString();
			for (int i = 0; i < power; i++) {
				r = r.times( Asigmas.get(symbol) );
			}
			sequence = sequence.substring(nextstreak.rawStringLength(), sequence.rawStringLength()); //+1 may be needed to get rid of comma
		}
		
		assert(Math.abs(r.times(ainf).get(0, 0)) == r.norm1() );
		return r.times(ainf).get(0, 0);
		
	}
	
		
}
