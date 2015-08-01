package hmm_sim;

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.PriorityQueue;
import java.util.Set;

import javax.swing.text.StyleContext.SmallAttributeSet;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModelMultipleObservations extends HankelSVDModelParent {
	
	//String formatting: 	"31:15,21:12", IDofSymbol:Count,IDofSymbol:Count ... etc
	
	private SymbolCounts fullData;
	private SymbolCounts prefixes;
	private SymbolCounts suffixes;
	
	private int numDimensions;
	private int basisSize;
	private SingularValueDecomposition svdOfH;
	
	public static void main(String[] args){
		HankelSVDModelMultipleObservations.test();
	}
		
	public static void test(){
		//String[] samples = {"1:2" , "2:1,1:1"};
		String[] samples = {"1:1","1:1,2:1","1:2,2:1","1:1,2:1,1:1", "2:1,1:2"};
		LinkedList<SymbolCountPair> l = new LinkedList<SymbolCountPair>();
		for (int i = 0; i < samples.length; i++) {
			String s = samples[i];
			SequenceOfSymbols seq = new SequenceOfSymbols(s);
			SymbolCountPair sc = new SymbolCountPair(1, seq);
			l.add(sc);
		}
		
		HankelSVDModelMultipleObservations h = new HankelSVDModelMultipleObservations(l, 15, 2);
	}
	
	public HankelSVDModelMultipleObservations(LinkedList<SymbolCountPair> scp, int basisSize, int numDimensions){
		this.numDimensions = numDimensions;
		this.basisSize = basisSize;
		this.fullData = new SymbolCounts(numDimensions);
		for (SymbolCountPair s: scp) {
			this.fullData.updateFrequency(s.getSymbol(), 1);
		}
		//System.out.println("FullData");
		//System.out.println(fullData);
		//System.out.println();
		
		SymbolCounts prefixes = getPrefixes(scp);
		SymbolCounts suffixes = getSuffixes(prefixes);
		
		Matrix Hlambda = this.buildHankelMultipleObservations(fullData, prefixes, suffixes, new SequenceOfSymbols( "") );
		
		this.svdOfH =  new SingularValueDecomposition(Hlambda);
	}
	
	private SymbolCounts getPrefixes(Iterable<SymbolCountPair> spp){

		SymbolCounts preFixes = new SymbolCounts(this.numDimensions);
		for (SymbolCountPair symbolCountPair : spp) {
			SequenceOfSymbols s = symbolCountPair.getSymbol();
			LinkedList<SequenceOfSymbols> prefixesFromSymbol = s.getPrefixesFromSequence();
			for (SequenceOfSymbols seq: prefixesFromSymbol) {
				preFixes.insertKeyTreatAsHashSet(seq); 
			}
		}
		
		SymbolCounts scReturn = new SymbolCounts(this.numDimensions);
		NavigableSet<SequenceOfSymbols> prefixesSorted = preFixes.incrKeySet();
		int i=0;
		for (SequenceOfSymbols s: prefixesSorted) {
			
			if (i >= basisSize){
				break;
			}
			else{
				scReturn.updateFrequency(s, preFixes.symbolToFrequency.get(s) );
				i++;
			}
		}
		/*
		System.out.println("Done choosing Prefixes: -HankelSVDModelMultipleObservations.getPrefixes");
		System.out.println("PreCut");
		System.out.println(preFixes);
		System.out.println("PostCut");
		System.out.println(scReturn);
		System.out.println();
		*/
		
		return scReturn;
		
	}

	private SymbolCounts getSuffixes(SymbolCounts prefixes){
		SymbolCounts suffixes = new SymbolCounts(this.numDimensions);
		for (SequenceOfSymbols prefix : prefixes.incrKeySet() ) {
			LinkedList<SequenceOfSymbols> L = prefix.getSuffixesOfSequence();
			for (SequenceOfSymbols sequenceOfSymbols : L) {
				suffixes.insertKeyTreatAsHashSet(sequenceOfSymbols); 	//suffix frequency is irrelevant in the way we choose the basis
			}
		}
		/*
		System.out.println("Done choosing Suffixes: -HankelSVDModelMultipleObservations.getSuffixes");
		System.out.println("Here they are: ");
		System.out.println(suffixes);
		System.out.println();
		*/
		
		return suffixes;
	}
	
	public Matrix buildHankelMultipleObservations(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes, SequenceOfSymbols X){
		double[][] hankel = new double[prefixes.incrKeySet().size()][suffixes.incrKeySet().size()];
		int[] freqCounter = determineTotalFrequencyIncludedPerMinKLength(dataCounts, prefixes, suffixes);
		/*System.out.println("Frequency Counter - buildHankelMultipleObsverations");
		System.out.println(Arrays.toString(freqCounter));
		System.out.println();
		*/
		int pC = 0;
		NavigableSet<SequenceOfSymbols> prefixKeySetSorted = prefixes.incrKeySet(); 
		NavigableSet<SequenceOfSymbols> suffixKeySetSorted = suffixes.incrKeySet(); 

		int L = getLengthOfLongestSequence(dataCounts);
		double[] sumOfSymbolLengthK = new double[L+1];
				
		for (SequenceOfSymbols prefkey: prefixKeySetSorted) {
			int pS = 0;
			for(SequenceOfSymbols suffkey: suffixKeySetSorted){
				SequenceOfSymbols t = SequenceOfSymbols.concatenateSymbols(prefkey, X);
				SequenceOfSymbols s = SequenceOfSymbols.concatenateSymbols(t, suffkey);
				
				if (dataCounts.getSymbolToFrequency().containsKey(s)){
					double tempDouble = (double) dataCounts.getSymbolToFrequency().get(s)/freqCounter[s.sequenceLength()];
					sumOfSymbolLengthK[s.sequenceLength()] += tempDouble;
					hankel[pC][pS] = tempDouble;
				}
				else{
					hankel[pC][pS] = 0;
				}
				
				pS++;
			}
			pC ++;
		}
		
		System.out.println("Verifying that F computes probabilities in the right way:");
		System.out.println( Arrays.toString(sumOfSymbolLengthK) );
		System.out.println();
		
		System.out.println("Printing out Hankel");
		Matrix H = new Matrix(hankel);
		H.print(5, 5);
		System.out.println("Prefixes");
		int[] prefArray = SequenceOfSymbols.getRawSequencesRowVector(prefixes);
		System.out.println("Number of prefixes: " + Integer.toString(prefArray.length));
		System.out.println(	Arrays.toString( prefArray )) ;
		System.out.println();
		System.out.println("Suffixes");
		int[] suffArray = SequenceOfSymbols.getRawSequencesRowVector(suffixes) ;
		System.out.println("Number of suffixes: " + Integer.toString(suffArray.length));
		System.out.println(	Arrays.toString( suffArray )) ;
		System.out.println();
		
		return H;
	}
	
	private int[] determineTotalFrequencyIncludedPerMinKLength(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes){
		int L = getLengthOfLongestSequence(dataCounts)+1;
		int[] freqCounter = new int[L];
		
		for (SequenceOfSymbols d: dataCounts.incrKeySet()) {
			for (int i = 0 ; i <= d.sequenceLength(); i++) {
				freqCounter[i] ++;
			}
		}
			
		/* Only including coverage of prefixes and suffixes. NA for borjas algorithm 
		HashSet<SequenceOfSymbols> seenSequences = new HashSet<SequenceOfSymbols>();
		for (SequenceOfSymbols prefkey: prefixes.incrKeySet() ){
			for(SequenceOfSymbols suffkey: suffixes.incrKeySet() ){
				SequenceOfSymbols s = SequenceOfSymbols.concatenate(prefkey, suffkey);
				if (seenSequences.contains(s) == false && dataCounts.getSymbolToFrequency().containsKey(s)){
					System.out.println(s);
					seenSequences.add(s);
					for (int i = 0 ; i <= s.sequenceLength(); i++) {
						freqCounter[i] ++;
					}
				}
			}
		}
		*/
		return freqCounter;
	}
	
	private int getLengthOfLongestSequence(SymbolCounts dataCounts){
		int l = 0;
		
		for (SequenceOfSymbols s: dataCounts.incrKeySet()) {
			int sl = s.sequenceLength();
			if (sl > l){
				l = sl;
			}
		}
		return l;
	}
	
	public QueryEngineMultipleObservations buildHankelBasedModelMultipleObservations(int base, int modelSize){
		SymbolCounts dataCounts = this.fullData;
		SymbolCounts prefixes = this.prefixes;
		SymbolCounts suffixes = this.suffixes;

		Matrix H = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, new SequenceOfSymbols(""));
				
		HashMap<String, Matrix> truncatedSVD = super.truncateSVD(H , modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
		
		HashMap<String, Matrix> XSigmas = new HashMap<String, Matrix>();
	
		for (int i = 1; i < this.numDimensions; i++) {
			String iString = Integer.toString(i);
			Matrix Hx = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, new SequenceOfSymbols( iString ) );
			Matrix Ax = pinv.times(Hx).times(sinv);
			XSigmas.put(iString, Ax);
		}
		
		double[][] h_L = new double[1][basisSize];
		int i = 0;
		for (SequenceOfSymbols s: this.prefixes.incrKeySet()) {
			int counts = this.prefixes.getSymbolToFrequency().get(s);
			h_L[0][i] = counts;
			i++;
		}
		i=0;
		
		int j = 0; 
		double[][] h_t = new double[basisSize][1];
		for (SequenceOfSymbols s: this.suffixes.incrKeySet()) {
			int counts = this.suffixes.getSymbolToFrequency().get(s);
			h_t[j][0] = counts;
			j++;
		}
		j=0;
		
		Matrix h_LS = new Matrix( h_L ).transpose();
		Matrix h_LT = new Matrix( h_t ).transpose();
				
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_LT);
	
		QueryEngineMultipleObservations q = new QueryEngineMultipleObservations(alpha_0, alpha_inf, XSigmas, base);
		
		//If Debugging Wanted
		//QueryEngine q = new QueryEngine(alpha_0, alpha_inf, Asigmas, maxExponent, base , pinv, sinv, truncatedSVD, this.svd);

		return q;
	}

	public synchronized void writeObject(java.io.ObjectOutputStream stream){}
	
	public void readObject(java.io.ObjectInputStream in){}

	@Override
	public QueryEngine buildHankelBasedModel(int base, int modelSize) {
		fill in
	}
}
