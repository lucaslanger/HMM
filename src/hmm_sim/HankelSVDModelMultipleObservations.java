package hmm_sim;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.PriorityQueue;

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
		String[] samples = {"1:2" , "2:1,1:1"};
		LinkedList<SymbolCountPair> l = new LinkedList<SymbolCountPair>();
		for (int i = 0; i < samples.length; i++) {
			String s = samples[i];
			SequenceOfSymbols seq = new SequenceOfSymbols(s);
			SymbolCountPair sc = new SymbolCountPair(1, seq);
			l.add(sc);
		}
		
		HankelSVDModelMultipleObservations h = new HankelSVDModelMultipleObservations(l, 3, 2);
	}
	
	public HankelSVDModelMultipleObservations(LinkedList<SymbolCountPair> scp, int basisSize, int numDimensions){
		this.numDimensions = numDimensions;
		this.basisSize = basisSize;
		this.fullData = new SymbolCounts(numDimensions);
		for (SymbolCountPair s: scp) {
			this.fullData.updateFrequency(s.getSymbol(), s.getCount());
		}
		SymbolCounts prefixes = getPrefixes(scp);
		SymbolCounts suffixes = getSuffixes(prefixes, scp);
		
		Matrix Hlambda = this.buildHankelMultipleObservations(fullData, prefixes, suffixes, "");
		
		this.svdOfH =  new SingularValueDecomposition(Hlambda);
	}
	
	private SymbolCounts getPrefixes(Iterable<SymbolCountPair> spp){

		SymbolCounts preFixes = new SymbolCounts(this.numDimensions);
		for (SymbolCountPair symbolCountPair : spp) {
			SequenceOfSymbols s = symbolCountPair.getSymbol();
			LinkedList<SequenceOfSymbols> prefixesFromSymbol = s.getPrefixesFromSequence();
			for (SequenceOfSymbols seq: prefixesFromSymbol) {
				preFixes.updateFrequency(seq, symbolCountPair.getCount() );
			}
		}
		
		SymbolCounts scReturn = new SymbolCounts(this.numDimensions);
		NavigableSet<SequenceOfSymbols> prefixesSorted = preFixes.descKeySet();
		int i=0;
		for (SequenceOfSymbols s: prefixesSorted) {
			if (i >= basisSize){
				break;
			}
			else{
				scReturn.updateFrequency(s, preFixes.getDataCount() );
				i++;
			}
		}
		
		return scReturn;
		
	}

	private SymbolCounts getSuffixes(SymbolCounts prefixes, Iterable<SymbolCountPair> spp){
		return null;
	}
	
	public Matrix buildHankelMultipleObservations(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes, String X){
		double[][] hankel = new double[prefixes.SortedKeys().size()][suffixes.SortedKeys().size()];
		int freqCounter = determineTotalFrequencyIncluded(dataCounts, prefixes, suffixes);
		
		int pC = 0;
		NavigableSet<SequenceOfSymbols> prefixKeySetSorted = prefixes.descKeySet(); 
		NavigableSet<SequenceOfSymbols> suffixKeySetSorted = suffixes.descKeySet(); 

		for (SequenceOfSymbols prefkey: prefixKeySetSorted) {
			int pS = 0;
			for(SequenceOfSymbols suffkey: suffixKeySetSorted){
				String s = prefkey + X + suffkey;
				if (dataCounts.getSymbolToFrequency().containsKey(s)){
					hankel[pC][pS] = (double) dataCounts.getSymbolToFrequency().get(s)/freqCounter;
				}
				else{
					hankel[pC][pS] = 0;
				}
				
				pS++;
			}
			pC ++;
		}
		Matrix H = new Matrix(hankel);
		return H;
	}
	
	public int determineTotalFrequencyIncluded(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes){
		int freqCounter = 0;
		for (SequenceOfSymbols prefkey: prefixes.descKeySet() ){
			for(SequenceOfSymbols suffkey: suffixes.descKeySet() ){
				SequenceOfSymbols s = SequenceOfSymbols.concatenate(prefkey, suffkey);
				freqCounter += dataCounts.getSymbolToFrequency().get(s);
			}
		}
		return freqCounter;
	}
	
	public QueryEngineMultipleObservations buildHankelBasedModelMultipleObservations(int base, int modelSize){
		SymbolCounts dataCounts = this.fullData;
		SymbolCounts prefixes = this.prefixes;
		SymbolCounts suffixes = this.suffixes;

		Matrix H = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, "");
				
		HashMap<String, Matrix> truncatedSVD = super.truncateSVD(H , modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
		
		HashMap<String, Matrix> XSigmas = new HashMap<String, Matrix>();
	
		for (int i = 1; i < this.numDimensions; i++) {
			String iString = Integer.toString(i);
			Matrix Hx = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, iString );
			Matrix Ax = pinv.times(Hx).times(sinv);
			XSigmas.put(iString, Ax);
		}
		
		double[][] h_L = new double[1][basisSize];
		int i = 0;
		for (SequenceOfSymbols s: this.prefixes.descKeySet()) {
			int counts = this.prefixes.getSymbolToFrequency().get(s);
			h_L[0][i] = counts;
			i++;
		}
		i=0;
		
		int j = 0; 
		double[][] h_t = new double[basisSize][1];
		for (SequenceOfSymbols s: this.suffixes.descKeySet()) {
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
