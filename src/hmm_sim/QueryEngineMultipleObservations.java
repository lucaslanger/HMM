package hmm_sim;


import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class QueryEngineMultipleObservations {
		
	private Matrix a0;
	private Matrix ainf;
	private HashMap<SequenceOfSymbols, Matrix> Asigmas;
	private int maxExponent;
	private int maxPower;
	
	private boolean debug;
	//Fields below are only useful for debugging
	
	private Matrix pinv;
	private Matrix sinv;
	private HashMap<String, Matrix> truncatedSVD;
	private SingularValueDecomposition originalSVD;
	
	private SymbolCounts prefixes;
	private SymbolCounts suffixes;
	private int base;
	
	public int getBase(){
		return this.base;
	}
	
	public int getModelSize(){
		return this.a0.getArrayCopy()[0].length;
	}
	
	public QueryEngineMultipleObservations(int maxExponent, Matrix a0, Matrix ainf, HashMap<SequenceOfSymbols, Matrix> Asigmas, int base, SymbolCounts prefixes, SymbolCounts suffixes){
		this.a0 = a0;
		this.ainf = ainf;
		this.Asigmas = Asigmas;
		this.prefixes = prefixes;
		this.suffixes = suffixes;
		this.base = base;
		
		this.maxExponent = maxExponent;
		this.maxPower = (int) Math.pow(base, maxExponent);
		
		this.debug = false;
	}
	
	public HashMap<String, Matrix> getTruncatedSVD() {
		return truncatedSVD;
	}

	public SymbolCounts getPrefixes() {
		return prefixes;
	}

	public SymbolCounts getSuffixes() {
		return suffixes;
	}

	public QueryEngineMultipleObservations(int maxExponent, Matrix a0, Matrix ainf, HashMap<SequenceOfSymbols, Matrix> Asigmas, int base, SymbolCounts prefixes, SymbolCounts suffixes, Matrix pinv, Matrix sinv, HashMap<String, Matrix> truncatedSVD, SingularValueDecomposition originalSVD){
		this.a0 = a0;
		this.ainf = ainf;
		this.Asigmas = Asigmas;
		this.prefixes = prefixes;
		this.suffixes = suffixes;
		this.base = base;
		
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

	public HashMap<SequenceOfSymbols, Matrix> getAsigmas() {
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

	public double probabilityQuery(SequenceOfSymbols sequence, boolean debug){
		Matrix r = this.a0;
		while(sequence.getSequence() != ""){
			
			SequenceOfSymbols nextstreak = sequence.getFirstStreak();
			int power = nextstreak.getStreakFromString();
			String symbol = nextstreak.getSymbolFromString();
		
			r = processQueryFixedSymbol(r, symbol, power, debug);
			
			if (sequence.rawStringLength() > nextstreak.rawStringLength()){
				sequence = sequence.substring(nextstreak.rawStringLength()+1, sequence.rawStringLength());
			}
			else{
				sequence = new SequenceOfSymbols("");
			}
		}
		
		assert(Math.abs(r.times(ainf).get(0, 0)) == r.norm1() );
		return r.times(ainf).get(0, 0);
		
	}
	
	private Matrix processQueryFixedSymbol(Matrix r, String symbol, int power, boolean debug) {
		int currentLimitingPower = this.maxPower;
		SequenceOfSymbols currentSequence = new SequenceOfSymbols( symbol + ":" + Integer.toString(currentLimitingPower)  );

		//System.out.println(Asigmas.keySet());
		while (power > 0){	
			while(currentLimitingPower <= power){
				r = r.times( this.Asigmas.get(currentSequence) );
				power = power - currentLimitingPower;
			}
			currentLimitingPower /= base;
			currentSequence = new SequenceOfSymbols( symbol + ":" + Integer.toString(currentLimitingPower) );
		}
		
		return r;
	}

	public double evaluateModel(LabyrinthGraph L, HashSet<SequenceOfSymbols> stringsToQuery, boolean topErrors){
		double e = 0;
				
		PriorityQueue<SequenceErrorPair> pq = new PriorityQueue<SequenceErrorPair>();
		
		for (SequenceOfSymbols sequenceOfSymbols : stringsToQuery) {
			double probQ = this.probabilityQuery(sequenceOfSymbols, false);
			double realProb =  L.determineRealProbabilityOfSequenceDoubleLoop(sequenceOfSymbols);
			
			double error = Math.pow( probQ - realProb, 2);
			//double error = Math.abs( probQ - realProb);
			pq.add( new SequenceErrorPair(sequenceOfSymbols, -1*error));
			e += error;
		}
		
		if (topErrors){
			for (int i = 0; i < 10; i++) {
				SequenceErrorPair a = pq.remove();
				System.out.println( a );
				System.out.println("Real: " + L.determineRealProbabilityOfSequenceDoubleLoop(a.getSeq()));
				System.out.println("Computed: " + this.probabilityQuery(a.getSeq(), false));;
				System.out.println();
				
			}
		}
		e = Math.sqrt(e);
		return e;
	}
	
	public void verifyProbabilityQueryCorrectness(){
		System.out.println("Making sure that probability query's are behaving as they are supposed to.");
		String[] tests = new String[]{"1:15", "2:13", "2:12,1:13"};
		for (String string : tests) {
			SequenceOfSymbols seq = new SequenceOfSymbols(string);
			System.out.println(seq);
			this.probabilityQuery(seq, true);
			System.out.println();
		}
	}
	
		
}
