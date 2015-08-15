package hmm_sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.PriorityQueue;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class HankelSVDModelMultipleObservations extends HankelSVDModelParent {
	
	//String formatting: 	"31:15,21:12", IDofSymbol:Count,IDofSymbol:Count ... etc
	
	private SymbolCounts fullData;
	private SymbolCounts prefixes;
	private SymbolCounts suffixes;
	
	private int numDimensions;
	private int basisSize;
	private HashMap<String, Matrix> svdOfH;
	private int rank;
	
	public int getRank(){
		return this.rank;
	}
	
	public static void main(String[] args){
		//HankelSVDModelMultipleObservations.initialtest();
		HankelSVDModelMultipleObservations.doubleLoopTest();
	}
	
	public HashSet<SequenceOfSymbols> generatePowerBase(int base, int maxPower, int numDimensions){
		return null;
	}
		
	public static void doubleLoopTest(){
		String workingFolder = "keySearchPacMan/MultipleObservationPlots/";
		//FlowControl.createFolder(workingFolder);
		
		int repetitions = 2;
		int numberOfTrajectories = 1000;
		int amountOfData = 1000;
		
		int numDimensions = 2;
		int base = 2; 
		
		int loop1 = 27;
		int loop2 = 27;
		int desiredHankelSize = (loop1+loop2)*3;
		int basisSize = 35;
		String dataSetFolder = workingFolder + "DataSets"+ loop1 + ":" + loop2+ "/";
			
		int[] modelSizes = new int[]{5,10,15,20,25,30,35,40,45,50};
		int maxPower = 243;
		int maxExponent = (int) (Math.log(maxPower)/Math.log(base));
		
		//Leave commented if datasets are already there!
		//FlowControl.createFolder(dataSetFolder);
		//generateDataSet(repetitions, dataSetFolder, desiredHankelSize, numberOfTrajectories, loop1, loop2);
		//
		
		Matrix Eavg = null;
		double[][] xaxis = new double[maxExponent+1][modelSizes.length];
		LabyrinthGraph L = LabyrinthGraph.multipleObservationDoubleLoop(workingFolder, desiredHankelSize, numberOfTrajectories, loop1, loop2);
		
		for (int r = 0; r < repetitions; r++) {
			File f = new File(dataSetFolder + "DataSet" + ",Rep:" + (r+1));
			SequenceOfSymbols[] seqsRead = HankelSVDModelMultipleObservations.readDataSetFromFile(f, amountOfData);
		
			LinkedList<SymbolCountPair> l = new LinkedList<SymbolCountPair>();
			
			boolean printSequences = false;
			for (SequenceOfSymbols s : seqsRead) {
				if (printSequences){
					System.out.print("Length:" + s.sequenceLength() + ", " + s + "	");
				}
				SymbolCountPair sc = new SymbolCountPair(1, s);
				l.add(sc);
			}
			if (printSequences){
				System.out.println();
			}
			
			HankelSVDModelMultipleObservations h = new HankelSVDModelMultipleObservations(l, basisSize, numDimensions);
			
			HashSet<SequenceOfSymbols> stq = makeStringsToQuery(h.prefixes, h.suffixes);
			System.out.println("Number of strings used when computing fhat v.s f: " + stq.size() + "\n");
			
			double[][] errors = new double[maxExponent+1][modelSizes.length];	
			
			for (int i = 0; i <= maxExponent; i++) {
					
				QueryEngineMultipleObservations[] qs = makeEnginesFromSamples(h, seqsRead, numDimensions, basisSize, base, i, modelSizes);
				
				for (int j = 0; j < qs.length; j++) {
					quickQueryTest(qs[j]);
					boolean debugErrors = false;
					double e = qs[j].evaluateModel(L, stq, debugErrors);
					
					errors[i][j] = e;
					xaxis[i][j] = modelSizes[j];
					
					//qs[i].verifyProbabilityQueryCorrectness();
				}
			}
			
			Matrix E = new Matrix(errors);
			if (Eavg == null){
				h.plotSingularValues(workingFolder, loop1, loop2, amountOfData);	//Plot the first set of singular Values -- HACKY
				Eavg = E;
			}
			else{
				Eavg = Eavg.plus(E);
			}
		}
		Eavg = Eavg.times(1.0/repetitions);
		
		System.out.println("Rows: exponents, Columns: modelSizes");
		Eavg.print(5, 5);
		
		
		String title = "Wall Color Predictions";
		String internalComment = "Lighter Curves --> Less Base System";
		OutputData.outputData(workingFolder + "errorModelSizesBase" + "," + amountOfData + "," + + loop1 + ":" + loop2, "Model Size", "Error norm2()", xaxis, Eavg.getArrayCopy(), title, internalComment);
	
	}
	
	private static void generateDataSet(int repetitions, String workingFolder, int desiredHankelSize, int numberOfTrajectories, int loop1, int loop2) {
		for (int i = 0; i < repetitions; i++) {
			LabyrinthGraph t = LabyrinthGraph.multipleObservationDoubleLoop(workingFolder, desiredHankelSize, numberOfTrajectories, loop1, loop2);
			SequenceOfSymbols[] seqs = t.getData();
			File filename = new File(workingFolder + "DataSet" + ",Rep:" + (i+1) );
			HankelSVDModelMultipleObservations.writeDataSetToFile(filename , seqs);
		}
		
	}

	public static void writeDataSetToFile(File filename, SequenceOfSymbols[] seqs){
		try{
			ObjectOutputStream oop = new ObjectOutputStream( new FileOutputStream(filename) );
			for (SequenceOfSymbols sequenceOfSymbols : seqs) {
				oop.writeObject(sequenceOfSymbols);
			}
			oop.close();
		}
		catch(Exception e){
			e.printStackTrace();
		}
	}
	
	public static SequenceOfSymbols[] readDataSetFromFile(File filename, int amountOfData){
		try {
			ObjectInputStream ois = new ObjectInputStream( new FileInputStream(filename) );
			SequenceOfSymbols[] seqs = new SequenceOfSymbols[amountOfData];
			for (int i = 0; i < amountOfData; i++) {
				seqs[i] = (SequenceOfSymbols) ois.readObject();
			}
			ois.close();
			return seqs;
			
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	public void plotSingularValues(String workingFolder, int loop1, int loop2, int amountOfData){
		double[][] yaxis = new double[1][rank];
		double[][] xaxis = new double[1][rank];

		for (int i = 0; i < rank; i++) {
			yaxis[0][i] = this.svdOfH.get("S").get(i, i);
			xaxis[0][i] = i;
		}
		String title = "Singular Values of H";
		String internalComment = "";
		OutputData.outputData(workingFolder + "SingularValues,"+ amountOfData + "," + loop1 + ":" + loop2, "i'th Singular Value", "", xaxis, yaxis, title, internalComment);

		
	}
	
	public static HashSet<SequenceOfSymbols> makeStringsToQuery(SymbolCounts prefixes, SymbolCounts suffixes){
		HashSet<SequenceOfSymbols> stringsToQuery = new HashSet<SequenceOfSymbols>();
		for (SequenceOfSymbols p : prefixes.incrKeySet()) {
			for (SequenceOfSymbols s : suffixes.incrKeySet()) {
				SequenceOfSymbols c = SequenceOfSymbols.concatenateSymbols(p, s);
				stringsToQuery.add(c);
			}
		}
		return stringsToQuery;
	}
	
	public static void quickQueryTest(QueryEngineMultipleObservations q){
		
		String[] tests = new String[]{"1:1", "2:1", "1:1,2:1", "1:5", "1:9"};
		for (String string : tests) {
			SequenceOfSymbols seq = new SequenceOfSymbols(string);
			//System.out.println("Query");
			//System.out.println(seq);
			//System.out.println(q.probabilityQuery(seq, false));
			//System.out.println();
		}
		
	}
	
	public static void initialtest(){
		//String[] samples = {"1:2" , "2:1,1:1"};
		String[] samples = {"1:1","1:1,2:1","1:2,2:1","1:1,2:1,1:1", "2:1,1:2"};
		LinkedList<SymbolCountPair> l = new LinkedList<SymbolCountPair>();
		for (int i = 0; i < samples.length; i++) {
			String s = samples[i];
			SequenceOfSymbols seq = new SequenceOfSymbols(s);
			SymbolCountPair sc = new SymbolCountPair(1, seq);
			l.add(sc);
		}
		
		int numDimensions = 2;
		int basisSize = 15;
		
		HankelSVDModelMultipleObservations h = new HankelSVDModelMultipleObservations(l, basisSize, numDimensions);
		int base = 2;
		int maxExp = 0;
		int modelSize = 6; 
		QueryEngineMultipleObservations a = h.buildHankelBasedModelMultipleObservations(h.fullData, h.prefixes, h.suffixes, base, maxExp, modelSize);
		for (String sa : samples) {
			System.out.println(a.probabilityQuery( new SequenceOfSymbols(sa), false ));
		}
	}
	
	
	public static QueryEngineMultipleObservations[] makeEnginesFromSamples(HankelSVDModelMultipleObservations h, SequenceOfSymbols[] seqs, int numDimensions, int basisSize, int base, int maxExponent, int[] modelSizes){
		
		QueryEngineMultipleObservations[] qEs = new QueryEngineMultipleObservations[modelSizes.length];
		for (int i = 0; i <= maxExponent; i++) {
			
			int c=0;
			for (int mS: modelSizes) {
			
				if (h.getRank() < mS){
					mS = h.getRank();
					//System.out.println("Truncation downgraded to trueModel size, too big a model!");
					//System.out.println();
				}
				
				qEs[c] = h.buildHankelBasedModelMultipleObservations(h.fullData, h.prefixes, h.suffixes, base, i ,mS);
				c++;
			}
			
		}
		
		return qEs;
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
		
		SymbolCounts prefixes = getPrefixes(this.fullData, scp);
		SymbolCounts suffixes = getSuffixes(this.fullData, prefixes);
		this.prefixes = prefixes;
		this.suffixes = suffixes;
		
		Matrix Hlambda = this.buildHankelMultipleObservations(fullData, prefixes, suffixes, new SequenceOfSymbols(""), true );
		this.rank = Hlambda.rank();

		//Hlambda.print(5, 5);
		System.out.println("Rank");
		System.out.println(this.getRank());
		System.out.println("Singular values");
		this.svdOfH = super.takeSVD(Hlambda);
		System.out.println( Arrays.toString( turnStoArray(svdOfH.get("S")) ) );
		System.out.println();
		
	
	}
	
	public static double[] turnStoArray(Matrix singularValues){
		double[] s = new double[singularValues.rank()];
		for (int i = 0; i < singularValues.rank(); i++) {
			s[i] = singularValues.get(i, i);
		}
		return s;
	}
	
	private SymbolCounts getPrefixes(SymbolCounts fullData, Iterable<SymbolCountPair> spp){

		SymbolCounts preFixes = new SymbolCounts(this.numDimensions);
		for (SymbolCountPair symbolCountPair : spp) {
			SequenceOfSymbols s = symbolCountPair.getSequence();
			LinkedList<SequenceOfSymbols> prefixesFromSymbol = s.getPrefixesFromSequence();
			for (SequenceOfSymbols seq: prefixesFromSymbol) {
				preFixes.updateFrequency(seq, 1);  
			}
		}
		
		PriorityQueue<SymbolCountPair> pq = new PriorityQueue<SymbolCountPair>();
		for (SequenceOfSymbols sequenceOfSymbols : preFixes.incrKeySet()) {
			SymbolCountPair scp = new SymbolCountPair(-1*preFixes.getSymbolToFrequency().get(sequenceOfSymbols), sequenceOfSymbols);
			pq.add( scp );
		}
		
		SymbolCounts scReturnPrefixes = new SymbolCounts(this.numDimensions);
		
		System.out.println();
		System.out.println("Prefix ordering info: ");
		int i=0;
		while(true){
			
			if (i >= basisSize){
				break;
			}
			else{
				SymbolCountPair scp = pq.remove();
				scReturnPrefixes.updateFrequency(scp.getSequence(), scp.getCount() );
				System.out.print(scp.getSequence() + " ");
				System.out.print(scp.getCount()*-1 + ", ");
				i++;
			}
		}
		System.out.println();
		System.out.println();
		
		scReturnPrefixes = verifyPrefixCompleteness(fullData, scReturnPrefixes);

		return scReturnPrefixes;
		
	}
	
	private SymbolCounts verifyPrefixCompleteness(SymbolCounts fullData, SymbolCounts scReturn) {
		int count = 0;
		SymbolCounts returnCounts = new SymbolCounts(scReturn.getNumDimensions());
		for (SequenceOfSymbols sp: scReturn.getSymbolToFrequency().keySet()) {
			for (SequenceOfSymbols s : sp.getPrefixesFromSequence()) {
				if( returnCounts.getSymbolToFrequency().containsKey( s ) == false){
					/*System.out.println("Prefix incomplete");
					System.out.println("In: " + sp.toString());
					System.out.println("Missing: " + s.toString());
					System.out.println();*/
					returnCounts.updateFrequency(s, fullData.getSymbolToFrequency().get(s) );
					count++;
					//Double check that sp itself gets included
				}
			}
		}
		System.out.println("Number of extra prefixes added: " + Integer.toString(count - scReturn.getSymbolToFrequency().size()) );
		System.out.println();
		
		return returnCounts;
	}
	
	private SymbolCounts getSuffixes(SymbolCounts trueData, SymbolCounts prefixes){
		SymbolCounts suffixes = new SymbolCounts(this.numDimensions);
		for (SequenceOfSymbols prefix : prefixes.incrKeySet() ) {
			LinkedList<SequenceOfSymbols> L = prefix.getSuffixesOfSequence();
			for (SequenceOfSymbols sequenceOfSymbols : L) {
				if (suffixes.getSymbolToFrequency().containsKey(sequenceOfSymbols) == false){
					if (trueData.getSymbolToFrequency().containsKey(sequenceOfSymbols) ){
						suffixes.updateFrequency(sequenceOfSymbols, trueData.getSymbolToFrequency().get(sequenceOfSymbols) ); 	//suffix frequency is irrelevant in the way we choose the basis
					}
					else{
						suffixes.updateFrequency(sequenceOfSymbols, 0);
					}
				}
			}
		}
		
		return suffixes;
	}
	
	private SequenceOfSymbols[] sortCounts(SymbolCounts counts){
		SequenceOfSymbols[] sortedCounts = new SequenceOfSymbols[counts.getSymbolToFrequency().size()];
		PriorityQueue<SymbolCountPair> pq = new PriorityQueue<SymbolCountPair>();
		for (SequenceOfSymbols seq : counts.getSymbolToFrequency().keySet()) {
			SymbolCountPair scp = new SymbolCountPair(counts.getSymbolToFrequency().get(seq)*-1, seq);
			pq.add(scp);
		}
		int i=0;
		
		while (pq.isEmpty() == false) {
			sortedCounts[i] = pq.remove().getSequence();
		
			i++;
		}

		return sortedCounts;
	}
	
	public Matrix buildHankelMultipleObservations(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes, SequenceOfSymbols X, boolean printOuts){
		double[][] hankel = new double[prefixes.incrKeySet().size()][suffixes.incrKeySet().size()];
		
		int pC = 0;
		SequenceOfSymbols[] prefixKeySetSorted = sortCounts(prefixes); 
	
		SequenceOfSymbols[] suffixKeySetSorted = sortCounts(suffixes); 

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
					/*if (printOuts){
						System.out.println("x: " + pC + " y: " + pS + " probability " + hankel[pC][pS] + " symbol: " + s.getRawSequence());
					}*/	
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
			System.out.println("Verifying that F computes probabilities in the right way. \n");
			System.out.println("Occurances: (How many strings of length >= K)");
			System.out.println( Arrays.toString(freqCounter) );
			System.out.println();
			System.out.println("Coverage: (How many strings of length >= K included -Ratio)");
			String[] t = new String[sumOfSymbolLengthK.length];
			for (int i = 0; i < sumOfSymbolLengthK.length; i++) {
				t[i] = Double.toString(sumOfSymbolLengthK[i]).substring(0, Math.min(Double.toString(sumOfSymbolLengthK[i]).length(), 4));
			}
			System.out.println( Arrays.toString(t) );
			System.out.println();
			
			System.out.println("Printing out Hankel");
			
			System.out.println("Number of prefixes: " + Integer.toString(prefixes.getSymbolToFrequency().keySet().size()));
			System.out.println("Number of suffixes: " + Integer.toString(suffixes.getSymbolToFrequency().keySet().size()));
			System.out.println();
			try{
				System.out.println("Prefixes");
				int[] prefArray = SequenceOfSymbols.getRawSequencesRowVector(prefixKeySetSorted);
				System.out.println(	Arrays.toString( prefArray )) ;
				System.out.println();
				System.out.println("Suffixes");
				int[] suffArray = SequenceOfSymbols.getRawSequencesRowVector(suffixKeySetSorted) ;
				System.out.println(	Arrays.toString( suffArray )) ;
				System.out.println();
			}
			catch(NumberFormatException e){
				System.out.println("Streaks too long to be printed out");
				System.out.println();
			}
			
		}
		Matrix H = new Matrix(hankel);
		
		/*System.out.println("Hankel Matrix:");
		System.out.println("Subscript " + X.getSequence());
		H.print(5, 5);	*/
		return H;
	}
	
	private int[] determineTotalFrequencyIncludedPerMinKLength(SymbolCounts dataCounts){
		int L = getLengthOfLongestSequence(dataCounts)+1;
		int[] freqCounter = new int[L];
		
		for (SequenceOfSymbols d: dataCounts.incrKeySet()) {
			freqCounter[d.sequenceLength()] += dataCounts.getSymbolToFrequency().get(d);
		}
			
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
	
	public QueryEngineMultipleObservations buildHankelBasedModelMultipleObservations(SymbolCounts dataCounts, SymbolCounts prefixes, SymbolCounts suffixes, int base, int maxExponent, int modelSize){

		//System.out.println("Chosen model Size: " + modelSize + "\n");

		Matrix H = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, new SequenceOfSymbols(""), false);
				
		HashMap<String, Matrix> truncatedSVD = super.truncateSVD(H , modelSize);
		
		Matrix di = pseudoInvDiagonal(truncatedSVD.get("S"));
		Matrix pinv = di.times(truncatedSVD.get("U").transpose());
		Matrix sinv = (truncatedSVD.get("VT")).transpose();
		
		/*
		System.out.println("Double check this computation!");

		printMatrixDimensions(H, "hankel");
		printMatrixDimensions(pinv, "pinv");
		printMatrixDimensions(sinv, "sinv");
		*/
		
		HashMap<SequenceOfSymbols, Matrix> XSigmas = new HashMap<SequenceOfSymbols, Matrix>();
	
		for (int i = 1; i <= this.numDimensions; i++) {
			for (int j = 0; j <= maxExponent; j++) {	
				int pow = (int) Math.pow(base, j);
				String power = Integer.toString(pow);
				String iString = Integer.toString(i) + ":" + power;
				SequenceOfSymbols seq = new SequenceOfSymbols(iString);
				
				Matrix Hx = buildHankelMultipleObservations(dataCounts, prefixes, suffixes, seq, false );
				Matrix Ax = pinv.times(Hx).times(sinv);
				XSigmas.put(seq, Ax);
			}
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
	
		
		//System.out.println(XSigmas.keySet());
		QueryEngineMultipleObservations q = new QueryEngineMultipleObservations(maxExponent, alpha_0, alpha_inf, XSigmas, base, prefixes, suffixes);
		
		//If Debugging Wanted
		//QueryEngine q = new QueryEngine(alpha_0, alpha_inf, Asigmas, maxExponent, base , pinv, sinv, truncatedSVD, this.svd);
		
		return q;
	}
	
	private Matrix extendMatrixWithZeroes(Matrix m, int extendTo) {
		assert extendTo >= m.getArrayCopy().length;
		
		double[][] n = new double[extendTo][m.getArrayCopy()[0].length];
		for (int i = 0; i < n.length; i++) {
			if (i < m.getArrayCopy().length){
				for (int j = 0; j < n[i].length; j++) {
					n[i][j] = m.get(i, j);
				}
			}
			else{
				for (int j = 0; j < n[i].length; j++) {
					n[i][j] = 0;
				}
			}
		}
		Matrix N = new Matrix(n);
		return N;
	}

	private static void printMatrixDimensions(Matrix h, String id){
		//System.out.println("Number of columns: " + h.getArrayCopy()[0].length);
		//System.out.println("Number of rows: " + h.getArrayCopy().length);
		System.out.println( "Dimensions of: " + id + " " + h.getArrayCopy().length + "x" + h.getArrayCopy()[0].length);
		System.out.println();
	}
	
	

	public synchronized void writeObject(java.io.ObjectOutputStream stream){}
	
	public void readObject(java.io.ObjectInputStream in){}

	@Override
	public QueryEngine buildHankelBasedModel(int base, int modelSize) {
		// TODO Auto-generated method stub
		return null;
	}

	
	
}
