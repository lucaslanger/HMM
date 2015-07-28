package hmm_sim;


import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.PriorityQueue;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModel implements Serializable{
	
	private double[] probabilities;
	private int basisSize;

	private SymbolInfo fullData;
	
	private SingularValueDecomposition svd;

	
	public static void main(String[] args){
	}

	
	public double[] getProbabilities() {
		return probabilities;
	}

	public SingularValueDecomposition getSvd() {
		return svd;
	}
	
	
	public int getBasisSize() {
		return basisSize;
	}
	
	public HankelSVDModel(){
		
	}

	public HankelSVDModel(SymbolInfo si, int basisSize){
		this.fullData = fullData;
		this.takeSVDMultipleObservations();
	}


	public HankelSVDModel(double[] probabilities , int basisSize){
		this.probabilities = probabilities;
		this.basisSize = basisSize;
		takeSVD();
	}
	
	public HankelSVDModel(double[] probabilities , int basisSize, SingularValueDecomposition s){
		this.probabilities = probabilities;
		this.basisSize = basisSize;
		this.svd = s;
	}
	

	private SymbolProbabilityPair[] getBasisMultipleObservations(){
		// Maybe improvement of the base is only due to how your picking basis
		PriorityQueue<SymbolProbabilityPair> pq = new PriorityQueue<SymbolProbabilityPair>();
		for( String s: this.fullData.getSymbolToProbability().keySet()){
			SymbolProbabilityPair spp = new SymbolProbabilityPair(fullData.getSymbolToProbability().get(s), s);
			pq.add(spp);
		}
		
		SymbolProbabilityPair[] d = new SymbolProbabilityPair[this.basisSize*2];
		for (int i = 0; i < basisSize*2; i++) {
			d[i] = pq.remove();
		}
		return d;
		
	}
	
	private void takeSVDMultipleObservations() {
		SymbolProbabilityPair[] spp = this.getBasisMultipleObservations();
		
		SymbolCountPair[] scp;
		try{
			HashMap<String, Integer> prefixes = new HashMap<String, Integer>();
			HashMap<String, Integer> suffixes = new HashMap<String, Integer>();
			
			for (int i = 0; i < spp.length; i++) {
				SymbolProbabilityPair sp = spp[i];
				for( String p: sp.getPrefixes()){
					instantiateOrIncrementHashMapCounter(prefixes, p, spp[i].getCount());
				}
			}
			
			
			SingularValueDecomposition svd = H.svd();
			this.svd = svd;
		
		}
		catch(Exception e){
			System.out.println("Problem making H in HankelSVDModel");
			e.printStackTrace();
		}
		
	}
	
	public Matrix buildHankelMultipleObservations(HashMap<String, Integer> dataCounts, HashMap<String, Integer> prefixes, HashMap<String, Integer> suffixes){
		double[][] hankel = new double[prefixes.keySet().size()][suffixes.keySet().size()];
		int freqCounter = determineTotalFrequencyIncluded(dataCounts, prefixes, suffixes);
		
		int pC = 0;
		for (String prefkey: prefixes.keySet() ){
			int pS = 0;
			for(String suffkey: suffixes.keySet()){
				String s = prefkey + suffkey;
				hankel[pC][pS] = (double) dataCounts.get(s)/freqCounter;
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
	
	public HashMap<String, Integer> instantiateOrIncrementHashMapCounter(HashMap<String,Integer> h, String k ,int i){
		if (h.containsKey(k)){
			int t = h.get(k);
			h.put(k, t + i);
		}
		else{
			h.put(k, i);
		}
		return h;
	}
	
	public QueryEngine buildHankelBasedModelMultipleObservations(int basisSize, int base, int modelSize){
		
		HashMap<String, Matrix> truncatedSVD = this.truncateSVD(modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
		
		HashMap<String, Matrix> XSigmas = new HashMap<String, Matrix>();
	
		for(String symbol: this.symbols){
			Matrix hx = buildHMultipleObservations(symbol);
			XSigmas.put(symbol, hx);
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < basisSize; i++) {
			h_L[i][0] = this.basis.getSymbolToIndex()[i];
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
	
	public void takeSVD(){
		try{
			Matrix H = buildH(0, basisSize);
			SingularValueDecomposition svd = H.svd();
			this.svd = svd;
		}
		catch(Exception e){
			e.printStackTrace();
		}
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
	
	public static double[] getDiagonalArray(Matrix m){
		double[] r = new double[m.getArrayCopy()[0].length];
		for (int i = 0; i < r.length; i++) {
			r[i] = m.get(i, i);
		}
		return r;
	}
	
	public Matrix buildH(int startingIndex, int endingIndex) throws Exception{		
		int hSize = (endingIndex - startingIndex);
		if ( (hSize + startingIndex)*2 > this.probabilities.length){
			throw new Exception("You asked for too large a Hankel Matrix. Increase the max durations recorded in your Environment.");
		} 
		else{
			double[][] hankel = new double[hSize][hSize];
			for (int i = 0; i < hSize; i++) {
				for (int j = 0; j < hSize; j++) {
					hankel[i][j] = this.probabilities[i+j+startingIndex];
				}
			}
			Matrix H = new Matrix(hankel);
			return H;
		}
		
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
	
	public int getRank(){
		return this.svd.getS().rank();
	}
	
	public Matrix getHankel(){
		return this.svd.getU().times(this.svd.getS()).times(this.svd.getV().transpose());
	}
	
	private synchronized void writeObject(java.io.ObjectOutputStream stream) throws java.io.IOException{
		stream.writeInt(this.probabilities.length);
		for (int i=0; i<this.probabilities.length; i++){
			stream.writeObject(this.probabilities[i]);
		}
		stream.writeInt(this.basisSize);
		stream.writeObject(this.svd);
	}
	
	private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException{
		int l = in.readInt();
		double[] probabilities = new double[l];
		for (int i = 0; i < l; i++) {
			probabilities[i] = (double) in.readObject();
		}
		this.probabilities = probabilities;
		this.basisSize = in.readInt();
		this.svd = (SingularValueDecomposition) in.readObject();
	}
	
	public static void testHankel(){
		Matrix Hbar = new Matrix( new double[][]{ {0,0.2,0.14}, {0.2,0.22,0.15}, {0.14,0.45,0.31} }).transpose();
		
		Matrix Ha = new Matrix(new double[][]{ {0.2,0.22,0.15},{0.22,0.19,0.13},{0.45,0.45,0.32} }).transpose();
		Matrix Hb = new Matrix(new double[][]{ {0.14,0.45,0.31}, {0.15,0.29,0.13}, {0.31,0.85,0.58} } ).transpose();
		Matrix hls = new Matrix(new double[][]{ {0, 0.2, 0.14} } );
		Matrix hpl = new Matrix(new double[][]{ {0, 0.2, 0.14} } ).transpose();
			
		Hbar.print(5,5);
		
		SingularValueDecomposition svd = Hbar.svd();
		Matrix p = svd.getU().times( svd.getS() );
		Matrix s = svd.getV().transpose();
		
		Matrix pinv = p.inverse();
		
		System.out.println("INVERSE TEST");
		p.times(pinv).print(5, 5);
		
		Matrix sinv = s.inverse();
		s.times(sinv).print(5, 5);
		
		Matrix Aa = pinv.times(Ha).times(sinv); 
		Matrix Ab = pinv.times(Hb).times(sinv); 
		
		Matrix alpha0 = hls.times(sinv);	//alpha0 row
		Matrix alphainf = pinv.times(hpl);	//alphainf column
		
		Matrix test1 = alpha0.times(Aa).times(alphainf);
		Matrix test2 = alpha0.times(Ab).times(alphainf);
		Matrix test3 = alpha0.times(Aa).times(Ab).times(alphainf);
		Matrix test4 = alpha0.times(Ab).times(Aa).times(alphainf);
		Matrix test5 = alpha0.times(Aa).times(Aa).times(alphainf);
		Matrix test6 = alpha0.times(Ab).times(Ab).times(alphainf);
		
		test1.print(5,5);
		test2.print(5,5);
		test3.print(5,5);
		test4.print(5,5);
		test5.print(5,5);
		test6.print(5,5);
	}
	
	public static Matrix pseudoInvDiagonal(Matrix m){
		double[][] a = m.getArrayCopy();
		for (int i = 0; i < a.length; i++) {
			if (a[i][i] != 0){
				a[i][i] = 1/a[i][i];
			}
			else{
				a[i][i] = 0;
			}
		}
		return new Matrix(a);
	}
		

}
