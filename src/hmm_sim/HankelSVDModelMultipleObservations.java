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
		h.buildHankelBasedModelMultipleObservations(h.fullData, h.prefixes, h.suffixes, 2, 5);
	}
	
	public HankelSVDModelMultipleObservations(LinkedList<SymbolCountPair> scp, int basisSize, int numDimensions){
		this.numDimensions = numDimensions;
		this.basisSize = basisSize;
		this.fullData = new SymbolCounts(numDimensions);
		for (SymbolCountPair s: scp) {
			for (SequenceOfSymbols seq : s.getSequence().getPrefixesFromSequence()) {
				this.fullData.updateFrequency(seq, 1);									//Count all prefixes as having appeared
			}
		}

		
		SymbolCounts prefixes = getPrefixes(scp);
		SymbolCounts suffixes = getSuffixes(prefixes);
		this.prefixes = prefixes;
		this.suffixes = suffixes;
		
		/*
		System.out.println("FullData");
		System.out.println(fullData);
		System.out.println();
		
		System.out.println("Prefixes");
		System.out.println(prefixes.incrKeySet());
		System.out.println("Suffixes");
		System.out.println(suffixes.incrKeySet());
		System.out.println();
		 */
		
		Matrix Hlambda = this.buildHankelMultipleObservations(fullData, prefixes, suffixes, new SequenceOfSymbols( ""), true );
		
		this.svdOfH =  new SingularValueDecomposition(Hlambda);
	}
	
	private SymbolCounts getPrefixes(Iterable<SymbolCountPair> spp){

		SymbolCounts preFixes = new SymbolCounts(this.numDimensions);
		for (SymbolCountPair symbolCountPair : spp) {
			SequenceOfSymbols s = symbolCountPair.getSequence();
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
		System.out.println(suffixes.incrKeySet());
		System.out.println();
		*/
		
		return suffixes;
	}
	
	public Matrix buildHankelMultipleObservations(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes, SequenceOfSymbols X, boolean printOuts){
		double[][] hankel = new double[prefixes.incrKeySet().size()][suffixes.incrKeySet().size()];
		
		int pC = 0;
		NavigableSet<SequenceOfSymbols> prefixKeySetSorted = prefixes.incrKeySet(); 
		NavigableSet<SequenceOfSymbols> suffixKeySetSorted = suffixes.incrKeySet(); 

		int[] freqCounter = determineTotalFrequencyIncludedPerMinKLength(dataCounts);
		
		int L = getLengthOfLongestSequence(dataCounts);
		double[] sumOfSymbolLengthK = new double[L+1];
		HashSet<String> seenSequences = new HashSet<String>();	
		
		for (SequenceOfSymbols prefkey: prefixKeySetSorted) {
			int pS = 0;
			for(SequenceOfSymbols suffkey: suffixKeySetSorted){
				SequenceOfSymbols t = SequenceOfSymbols.concatenateSymbols(prefkey, X);
				SequenceOfSymbols s = SequenceOfSymbols.concatenateSymbols(t, suffkey);
				
				if (dataCounts.getSymbolToFrequency().containsKey(s)){
					double occurancesOfS = (double) dataCounts.getSymbolToFrequency().get(s);
					double numStringsLengthSuffLarge = freqCounter[s.sequenceLength()];
					hankel[pC][pS] = occurancesOfS/numStringsLengthSuffLarge;
					
					if(printOuts && seenSequences.contains(s.getSequence()) == false ){
						//System.out.println(s);
						//System.out.println(occurancesOfS/numStringsLengthSuffLarge);
						//System.out.println();
					
						sumOfSymbolLengthK[s.sequenceLength()] += occurancesOfS/numStringsLengthSuffLarge;
						seenSequences.add(s.getSequence());
					}
				}
				else{
					hankel[pC][pS] = 0;
				}
				
				pS++;
			}
			pC ++;
		}
		if (printOuts){
			System.out.println("Verifying that F computes probabilities in the right way:");
			System.out.println( Arrays.toString(freqCounter) );
			System.out.println( Arrays.toString(sumOfSymbolLengthK) );
			System.out.println();
			
			System.out.println("Printing out Hankel");
		
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
		}
		Matrix H = new Matrix(hankel);
		H.print(5, 5);
		return H;
	}
	
	private int[] determineTotalFrequencyIncludedPerMinKLength(SymbolCounts dataCounts){
		int L = getLengthOfLongestSequence(dataCounts)+1;
		int[] freqCounter = new int[L];
		
		for (SequenceOfSymbols d: dataCounts.incrKeySet()) {
			freqCounter[d.sequenceLength()] += dataCounts.getSymbolToFrequency().get(d);
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
	
	public QueryEngineMultipleObservations buildHankelBasedModelMultipleObservations(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes, int base, int modelSize){

		Matrix H = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, new SequenceOfSymbols(""), false);
				
		HashMap<String, Matrix> truncatedSVD = super.truncateSVD(H , modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
		
		HashMap<String, Matrix> XSigmas = new HashMap<String, Matrix>();
	
		for (int i = 1; i < this.numDimensions; i++) {
			String iString = Integer.toString(i);
			Matrix Hx = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, new SequenceOfSymbols( iString ), false );
			Matrix Ax = pinv.times(Hx).times(sinv);
			XSigmas.put(iString, Ax);
		}
		
		int widthOfH = H.getArrayCopy()[0].length;
		int heightOfH = H.transpose().getArrayCopy()[0].length;
		
		double[][] h_W = new double[1][widthOfH];
		for (int j = 0; j < widthOfH; j++) {
			h_W[0][j] = H.get(0, j);
		}
		
		int j = 0; 
		double[][] h_H = new double[heightOfH][1];
		for (int i = 0; i < heightOfH; i++) {
			h_H[i][0] = H.get(i, 0);
		}
		
		Matrix H_W = new Matrix( h_W );
		Matrix H_H = new Matrix( h_H );
				
		Matrix alpha_0 = H_W.times(sinv);
		Matrix alpha_inf = pinv.times(H_H);
	
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
