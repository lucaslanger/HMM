package hmm_sim;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.PriorityQueue;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModelMultipleObservations extends HankelSVDModelParent {
	
	//String formatting: 	"31:15,21:12", IDofSymbol:Count,IDofSymbol:Count ... etc
	
	private SymbolCounts fullData;
	private HashMap<String, Integer> prefixes;
	private HashMap<String, Integer> suffixes;
	
	private int numDimensions;
	private int basisSize;
	private SingularValueDecomposition svdOfH;
	
	public static void main(String[] args){
		
	}
		
	public void test(){
		String[] samples = {"1:2" , "2:1,1:1"};
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
		this.fullData = new SymbolCounts(numDimensions);
		for (SymbolCountPair s: scp) {
			this.fullData.updateFrequency(s.getSymbol(), s.getCount());
		}
		HashMap<String, Integer> prefixes = getPrefixes(scp);
		HashMap<String, Integer> suffixes = getSuffixes(prefixes, scp);
		
		Matrix Hlambda = this.buildHankelMultipleObservations(fullData, prefixes, suffixes, "");
		
		this.svdOfH =  new SingularValueDecomposition(Hlambda);
	}
	
	private HashMap<String, Integer> getPrefixes(Iterable<SymbolCountPair> spp){
		return null;
	}
	
	private HashMap<String, Integer> getSuffixes(HashMap<String, Integer> prefixes, Iterable<SymbolCountPair> spp){
		return null;
	}
	
	public Matrix buildHankelMultipleObservations(SymbolCounts dataCounts, HashMap<String, Integer> prefixes, HashMap<String, Integer> suffixes, String X){
		double[][] hankel = new double[prefixes.keySet().size()][suffixes.keySet().size()];
		int freqCounter = determineTotalFrequencyIncluded(dataCounts, prefixes, suffixes);
		
		int pC = 0;
		for (String prefkey: prefixes.keySet() ){
			int pS = 0;
			for(String suffkey: suffixes.keySet()){
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
	
	public int determineTotalFrequencyIncluded(SymbolCounts dataCounts, HashMap<String, Integer> prefixes, HashMap<String, Integer> suffixes){
		int freqCounter = 0;
		for (String prefkey: prefixes.keySet() ){
			for(String suffkey: suffixes.keySet()){
				String s = prefkey + suffkey;
				freqCounter += dataCounts.getSymbolToFrequency().get(s);
			}
		}
		return freqCounter;
	}
	
	public QueryEngine buildHankelBasedModel(int base, int modelSize){
		SymbolCounts dataCounts = this.fullData;
		HashMap<String, Integer > prefixes = this.prefixes;
		HashMap<String, Integer > suffixes = this.prefixes;

		
		Matrix H = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, "");
				
		SingularValueDecomposition truncatedSVD = super.truncateSVD(H , modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.getS());
		Matrix pinv = di.times(truncatedSVD.getU().transpose());
		Matrix sinv = (truncatedSVD.getV().transpose()).transpose();
		
		HashMap<String, Matrix> XSigmas = new HashMap<String, Matrix>();
	
		for (int i = 1; i < this.numDimensions; i++) {
			String iString = Integer.toString(i);
			Matrix Hx = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, iString );
			Matrix Ax = pinv.times(Hx).times(sinv);
			XSigmas.put(iString, Ax);
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < prefixes.keySet().size(); i++) {
			hl[0][i] = prefixes.get(key);
		}
		
		Matrix h_LS = new Matrix( h_L ).transpose();
		Matrix h_PL = h_LS.transpose();
				
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_PL);
	
		QueryEngine q = new QueryEngine(alpha_0, alpha_inf, Xsigmas, maxExponent, base);
		
		//If Debugging Wanted
		//QueryEngine q = new QueryEngine(alpha_0, alpha_inf, Asigmas, maxExponent, base , pinv, sinv, truncatedSVD, this.svd);

		return q;
	}
	
	
	public synchronized void writeObject(java.io.ObjectOutputStream stream){}
	
	public void readObject(java.io.ObjectInputStream in){}
}
