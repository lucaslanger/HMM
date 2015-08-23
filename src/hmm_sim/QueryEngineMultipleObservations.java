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
	
	private boolean customBase = false;
	
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

	public QueryEngineMultipleObservations(Matrix alpha_0, Matrix alpha_inf,
			HashMap<SequenceOfSymbols, Matrix> xSigmas, SymbolCounts prefixes2,
			SymbolCounts suffixes2) {
		
		this.a0 = alpha_0;
		this.ainf = alpha_inf;
		this.Asigmas = xSigmas;
		this.prefixes = prefixes2;
		this.suffixes = suffixes2;
		this.customBase = true;
		
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
		if (this.customBase){
			return this.customProbabilityQuery(sequence);
		}
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
	
	private double customProbabilityQuery(SequenceOfSymbols sequence) {
		if (sequence.getSequence().equals("")){return a0.times(ainf).get(0, 0);}
		Matrix r = this.a0;
		HashSet<String> k = HeuristicsForPickingBase.seqHashSetToStringOfSymbols(Asigmas.keySet());
		try {
			StringIntPair opt = HeuristicsForPickingBase.computeOptimalCompositionsNeccesaryCoinStyle(k, sequence.getRawSequence(), "Min");
			String[] breakdown = opt.getS().split(",");
			for (String string : breakdown) {
				SequenceOfSymbols s = SequenceOfSymbols.fullStringToCompressed(string);
				r = r.times(this.Asigmas.get(s));
			}
			r = r.times(this.ainf);
	
			
			if (Math.abs(r.get(0, 0)) != r.norm1()){
				r.print(5, 5);
				System.out.println("Problem with matrix multiplication when querying!");
			}
			return r.get(0, 0);
		
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("issue computing optimal ordering! must fix!");
			return -100000;
		}
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
