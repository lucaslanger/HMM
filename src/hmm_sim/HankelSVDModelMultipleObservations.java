package hmm_sim;

import java.util.HashMap;
import java.util.LinkedList;
import java.util.PriorityQueue;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModelMultipleObservations {
	
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
	
	public QueryEngine buildHankelBasedModel(int basisSize, int base, int modelSize){
		
		HashMap<String, Matrix> truncatedSVD = this.truncateSVD(modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
		
		/*System.out.println("Testing inverses");
		Matrix test1 = pinv.times(truncatedSVD.get("U").times(truncatedSVD.get("S")));
		Matrix test2 = truncatedSVD.get("VT").times(sinv);
		System.out.println(Arrays.toString( getDiagonalArray(test1) ));
		System.out.println(Arrays.toString( getDiagonalArray(test2) ));
		*/
		
		int maxExponent = (int) Math.floor((Math.log( (this.probabilities.length/2) - basisSize)/Math.log(base))) ; 

		Matrix[] H_Matrices = new Matrix[maxExponent+1];
				
		int freq;
		Matrix h;
		//System.out.println("Building queryEngine");
		for (int l = 0; l <= maxExponent; l++) {
			freq = (int) Math.pow(base,l);
			try {
				h = this.buildH(freq, freq+basisSize);
				H_Matrices[l] = h;
			} catch (Exception e) {
				System.out.println("Problem Building Model when creating Hankel");
				e.printStackTrace();
				return null;
			}
		}
		
		Matrix Asigmas[] = new Matrix[maxExponent+1];
		Matrix t;

		for (int i = 0; i <= maxExponent; i++) {
			t = pinv.times(H_Matrices[i]).times( sinv );
			Asigmas[i] = pinv.times(H_Matrices[i]).times( sinv );
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < basisSize; i++) {
			h_L[i][0] = this.probabilities[i];
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
	
	public HashMap<String, Matrix> truncateSVD(int nStates){
		boolean debug = false;
		
	    Matrix U = this.svd.getU();
	    Matrix S = this.svd.getS();
	    Matrix V = this.svd.getV();
	    
	    double[][] utemp = U.getArrayCopy();
	    double[][] utrunc = new double[utemp.length][nStates];
	    for (int i = 0; i < utrunc.length; i++) {
	    	for (int j = 0; j < nStates; j++) {
	    		utrunc[i][j] = utemp[i][j];
			}
			
		}
	    Matrix Utrunc = new Matrix(utrunc);
	    
	    double[][] stemp = S.getArrayCopy();
	    double[][] strunc = new double[nStates][nStates];
	    for (int i = 0; i < nStates; i++) {
	    	for (int j = 0; j < nStates; j++) {
	    		strunc[i][j] = stemp[i][j];
			}
		}
	    Matrix Strunc = new Matrix(strunc);
	    
	    double[][] vtemp = V.getArrayCopy();			//Double check to make sure this isnt going wrong
	    double[][] vtrunc = new double[utemp.length][nStates];
	    for (int i = 0; i < vtrunc.length; i++) {
	    	for (int j = 0; j < nStates; j++) {
	    		vtrunc[i][j] = vtemp[i][j];
			}
		}
	    Matrix Vtrunc = new Matrix(vtrunc).transpose();

	    HashMap<String, Matrix> r = new HashMap<String, Matrix>();
	    r.put("U", Utrunc);
	    r.put("VT", Vtrunc);
	    r.put("S", Strunc);
	    
	    if (debug) {
	    	System.out.println("Before trunc");
	    	System.out.print("Size");
	    	System.out.println(S.getArrayCopy()[0].length);
	   	    //U.print(5, 5);
	   	    S.print(5, 5);
	   	    
	   	    //V.transpose().print(5, 5);
	   	    
	   	    System.out.println("After trunc");
	   	    System.out.print("Size");
	    	System.out.println(Strunc.getArrayCopy()[0].length);
	   	    //Utrunc.print(5, 5);
	   	    Strunc.print(5, 5);
	   	    //Vtrunc.print(5, 5);
		}    
	    
	    return r;
	    
	}
}
