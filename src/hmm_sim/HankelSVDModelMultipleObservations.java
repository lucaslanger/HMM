package hmm_sim;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.PriorityQueue;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModelMultipleObservations extends HankelSVDModelParent {
	
	private SymbolCounts fullData;
	private int numDimensions;
	private int basisSize;
	private SingularValueDecomposition svdOfH;
	
	public static void main(String[] args){
		String[] samples = {"aa" , "ba"};
		LinkedList<SymbolCountPair> l = new LinkedList<SymbolCountPair>();
		for (int i = 0; i < samples.length; i++) {
			String s = samples[i];
			SymbolCountPair sc = new SymbolCountPair(1, s);
			l.add(sc);
		}
		
		HankelSVDModelMultipleObservations h = new HankelSVDModelMultipleObservations(l, 3, 2);
	}
	
	public HankelSVDModelMultipleObservations(LinkedList<SymbolCountPair> scp, int basisSize, int numDimensions){
		this.numDimensions = numDimensions;
		this.basisSize = basisSize;
		fullData = new SymbolCounts(numDimensions);
		for (SymbolCountPair s: scp) {
			fullData.updateFrequency(s.getSymbol(), s.getCount());
		}
		svdOfH = this.takeSVDMultipleObservations();
	}
	
	
	private SymbolCountPair[] getBasisMultipleObservations(){
		// Maybe improvement of the base is only due to how your picking basis
		PriorityQueue<SymbolCountPair> pq = new PriorityQueue<SymbolCountPair>();
		for( String s: this.fullData.getSymbolToProbability().keySet()){
			SymbolCountPair spp = new SymbolCountPair(fullData.getSymbolToProbability().get(s), s);
			pq.add(spp);
		}
		
		SymbolCountPair[] d = new SymbolCountPair[this.basisSize*2];
		for (int i = 0; i < basisSize*2; i++) {
			d[i] = pq.remove();
		}
		return d;
		
	}
	
	private HashMap<String, Integer> getPrefixes(SymbolCountPair[] spp){
		return null;
	}
	
	private HashMap<String, Integer> getSuffixes(HashMap<String, Integer> prefixes, SymbolCountPair[] spp){
		return null;
	}
	
	public Matrix buildHankelMultipleObservations(HashMap<String, Integer> dataCounts, HashMap<String, Integer> prefixes, HashMap<String, Integer> suffixes, String X){
		double[][] hankel = new double[prefixes.keySet().size()][suffixes.keySet().size()];
		int freqCounter = determineTotalFrequencyIncluded(dataCounts, prefixes, suffixes);
		
		int pC = 0;
		for (String prefkey: prefixes.keySet() ){
			int pS = 0;
			for(String suffkey: suffixes.keySet()){
				String s = prefkey + X + suffkey;
				if (dataCounts.containsKey(s)){
					hankel[pC][pS] = (double) dataCounts.get(s)/freqCounter;
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
	
	public int determineTotalFrequencyIncluded(HashMap<String, Integer> dataCounts, HashMap<String, Integer> prefixes, HashMap<String, Integer> suffixes){
		int freqCounter = 0;
		for (String prefkey: prefixes.keySet() ){
			for(String suffkey: suffixes.keySet()){
				String s = prefkey + suffkey;
				freqCounter += dataCounts.get(s);
			}
		}
		return freqCounter;
	}
	
	public QueryEngine buildHankelBasedModelMultipleObservations(HashMap<String, Integer> dataCounts, HashMap<String, Integer> prefixes, HashMap<String, Integer> suffixes, int base, int modelSize){
		
		Matrix H = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, "");
				
		HashMap<String, Matrix> truncatedSVD = this.truncateSVDMultipleObservations(H , modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
		
		HashMap<String, Matrix> XSigmas = new HashMap<String, Matrix>();
	
		for(String symbol: this.symbols){
			Matrix Hx = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, symbol);
			Matrix Ax = pinv.times(Hx).times(sinv);
			XSigmas.put(symbol, Ax);
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < prefixes.keySet().size(); i++) {
			hl[0][i] = prefixes.get(key);
		}
		
		Matrix h_LS = new Matrix( h_L ).transpose();
		Matrix h_PL = h_LS.transpose();
				
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_PL);
	
		QueryEngine q = new QueryEngine(alpha_0, alpha_inf, Asigmas, maxExponent, base);
		
		//If Debugging Wanted
		//QueryEngine q = new QueryEngine(alpha_0, alpha_inf, Asigmas, maxExponent, base , pinv, sinv, truncatedSVD, this.svd);

		return q;
	}
	
	private HashMap<String, Matrix> truncateSVDMultipleObservations(Matrix h,
			int modelSize) {
		// TODO Auto-generated method stub
		return null;
	}
	
	
	public synchronized void writeObject(java.io.ObjectOutputStream stream){}
	
	public void readObject(java.io.ObjectInputStream in){}
}
