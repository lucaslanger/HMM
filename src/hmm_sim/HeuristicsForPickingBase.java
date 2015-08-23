package hmm_sim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.Set;

import javax.naming.BinaryRefAddr;

public class HeuristicsForPickingBase {
	
	public static void main(String[] args){
		//System.out.println("dsfsf".substring(1, 3));
		HeuristicsForPickingBase.testDPSelection();
	}
	public static void testDPSelection(){ 
		String[] t = new String[]{"1:1,2:2,1:1,2:1,1:1"};
			
			
		SequenceOfSymbols[] seqs = new SequenceOfSymbols[t.length];
		int i = 0;
		for(String s: t){
			seqs[i] = new SequenceOfSymbols(s);
			i++;
		}
		
		int numDimensions = 1;
		
		HashSet<String> base = new HashSet<String>();
		base.add("1");
		base.add("2");
		base.add("121");
		base.add("21");
				
		for (SequenceOfSymbols s: seqs) {
			try {
				StringIntPair opt = computeOptimalCompositionsNeccesaryCoinStyle(base, s.getRawSequence(), "Min" );
				System.out.println(s.getRawSequence());
				System.out.println(opt.getS());
				System.out.println(opt.getI());
				System.out.println();
			} catch (Exception e) {	e.printStackTrace();break;}
		}
		
	}
	
	public static void selectBaseTest(){
		String[] stemp = new String[]{"1:30", "1:30","1:60", "2:18"};
		SequenceOfSymbols[] seqs = new SequenceOfSymbols[stemp.length];
		int i=0;
		for (String s : stemp) {
			seqs[i] = new SequenceOfSymbols(s); 
			i++;
		}
		
		Map<SequenceOfSymbols, Integer> m = sequenceDataToCounts(seqs); 
		
		HashSet<SequenceOfSymbols> a = HeuristicsForPickingBase.chooseBaseFromData(m, 5, 2);

	}
	
	public static Map<SequenceOfSymbols, Integer> sequenceDataToCounts(SequenceOfSymbols[] seqs){
		SymbolCounts sc = new SymbolCounts();
		for (SequenceOfSymbols sequenceOfSymbols : seqs) {
			sc.updateFrequency(sequenceOfSymbols, 1);
		}
		return sc.symbolToFrequency;
	}
	

	public static StringIntPair computeOptimalCompositionsNeccesaryCoinStyle(HashSet<String> base, String s, String type){
		if (s.equals("")){
			return new StringIntPair("", 0);
		}
		
		HashMap<Integer, ArrayList<String>> possibleterminations = new HashMap<Integer, ArrayList<String>>();
		for (int i = 0; i < s.length(); i++) {
			for (String string : base) {
				if (string.equals("") == false){
					int stretchIndex = i + string.length() - 1 ;
					if ( stretchIndex < s.length() && string.equals(s.substring(i, stretchIndex + 1) )){
						if (possibleterminations.containsKey(stretchIndex)){
							ArrayList<String> t = possibleterminations.get(stretchIndex);
							t.add(string);
							possibleterminations.put(stretchIndex, t);
						}
						else{
							ArrayList<String> t = new ArrayList<String>();
							t.add(string);
							possibleterminations.put(stretchIndex, t);
						}
					}
				}
			}
			if (possibleterminations.containsKey(i) == false){
				System.out.println("Problem at index i");
				System.out.println(i);
				System.out.println(base);
				System.out.println(s);
			}
		}
		/*
		System.out.println("String: " + s);
		System.out.println();
		for (int i = 0; i < s.length(); i++) {
			System.out.println(i);
			for (String string : possibleterminations.get(i)) {
				System.out.println(string);
			}
			System.out.println();
		}
		*/
		
		
		if (type == "Min"){
			StringIntPair[] bestSolutions = new StringIntPair[s.length()];
			
			bestSolutions[0] = new StringIntPair( Character.toString(s.charAt(0)), 1);
			for (int i = 1; i < bestSolutions.length; i++) {
				StringIntPair best = null;
				boolean init = false;
				for (String seq : possibleterminations.get(i)) {
					if (i-seq.length() == -1){
						best = new StringIntPair(seq, 1);
						init = true;
					}
					else{
						int score = bestSolutions[i-seq.length()].getI() + 1;
						if (init == false || score < best.getI()){
							String sol = bestSolutions[i-seq.length()].getS() + "," + seq;
							best = new StringIntPair(sol, score);
							init = true;
						}
					}
				}
				bestSolutions[i] = best;
			}
			return bestSolutions[bestSolutions.length-1];
		}
		else{
			System.out.println("Type needs to be min right now!");
			return null;
		}
		
		
	}
	
	//Dynamic Programming algorithm first go, too slow
	public static StringIntPair computeOptimalCompositionsNeccesaryOldStyle(HashSet<String> base, String s, String type) throws Exception{
		if (s.equals("")){
			return new StringIntPair("", 0);
		}
		HashMap<Integer, ArrayList<String>> possiblePlugins = new HashMap<Integer, ArrayList<String>>();
		for (int i = 0; i < s.length(); i++) {
			for (String string : base) {
				if (string.equals("") == false){
					int stretchIndex = i + string.length() - 1 ;
					if ( stretchIndex < s.length() && string.equals(s.substring(i, stretchIndex + 1) )){
					
						if (possiblePlugins.containsKey(i)){
							ArrayList<String> t = possiblePlugins.get(i);
							t.add(string);
							possiblePlugins.put(i, t);
						}
						else{
							ArrayList<String> t = new ArrayList<String>();
							t.add(string);
							possiblePlugins.put(i, t);
						}
					}
				}
			}
			if (possiblePlugins.containsKey(i) == false){
				System.out.println("Problem at index i");
				System.out.println(i);
				System.out.println(base);
				System.out.println(s);
			}
		}
		int[][] compositionsNeccesary = new int[s.length()][s.length()];
		HashMap<String, String > bestCompositions = new HashMap<String, String>();
		
		return recursiveComputationsNeccessary(possiblePlugins, compositionsNeccesary, bestCompositions, base, s, 0, s.length()-1, type);
		
		
	}
	
	private static StringIntPair recursiveComputationsNeccessary(HashMap<Integer, ArrayList<String>> possiblePlugins, int[][] compositionsNeccesary, HashMap<String, String> bestCompositions, HashSet<String> base, String s, int i, int j, String type) throws Exception{
		
		if(i==j){
			String c = Character.toString(s.charAt(i));
			if(base.contains(  c ) ){
				String s1 = Integer.toString(i) + ":" +  j;
				bestCompositions.put(s1, c);
				return new StringIntPair(c,1);
			}
			else{
				System.out.println("Basic Symbol not included!" );
				System.out.println( s.charAt(i) );
				System.out.println(base);
				throw new Exception();
			}
		}
		
		else if(compositionsNeccesary[i][j] != 0 ){
			String s1 = Integer.toString(i) + ":" +  j;
			return new StringIntPair(bestCompositions.get(s1), compositionsNeccesary[i][j]);
		}
		else{
			int best = 0;
			String minBaseString = "";
			String bestComp = "";
			int bestH = -1;
			boolean b = false;
			for (int h = i; h <= j; h++) {
				for (String string : possiblePlugins.get(h)) {
					if (string.equals("")){throw new Exception();}
					
					if (string.length() - 1 + h > j){
						//donothing
					}else{
						StringIntPair cL;
						if(h != i){
							cL = recursiveComputationsNeccessary(possiblePlugins, compositionsNeccesary, bestCompositions, base, s, i, h-1, type);
						}else{
							cL = new StringIntPair("", 0);
						}
						//new addition to speed things up
						if (b != false && (cL.getI() + 1 >= best && type=="Min")){
							continue;
						}
						
						StringIntPair cR;
						if(h + string.length() <= j){
							cR = recursiveComputationsNeccessary(possiblePlugins, compositionsNeccesary, bestCompositions, base, s, h + string.length(), j, type);
						}
						else{
							cR = new StringIntPair("", 0);
						}
						
						int sum = cL.getI() + 1 + cR.getI();
						
						if (b == false || (sum < best && type=="Min") || (sum > best && type=="Max")){
							bestH = h;
							best = sum;
							minBaseString = string;
							bestComp = "";
							if (cL.getS().equals("") == false){
								bestComp = bestComp + cL.getS();
							} 
							bestComp = bestComp + minBaseString +  "," ;
							if (cR.getS().equals("") == false){
								bestComp = bestComp + cR.getS() ;
							} 
							b = true;
							
						}
					}
				
				}
			}
			if(b== false){
				throw new Exception();
			}
			
			String s1 = Integer.toString(i) +  ":" + j;
			compositionsNeccesary[i][j] = best;
			bestCompositions.put(s1, bestComp);
	
			return new StringIntPair(bestComp, best);
		}
	}
	
	
	//WORKING ON SUBSTRING STUFF
	
	public static HashSet<SequenceOfSymbols> chooseBaseFromData(Map<SequenceOfSymbols, Integer> m, int maxBaseSize, int numDimensions){
		HashSet<String> currentBase = new HashSet<String>();
		HashMap<SequenceOfSymbols, Integer> currentBestDecomposition = new HashMap<SequenceOfSymbols, Integer>();
		
		for (int i = 1; i <= numDimensions; i++) {
			currentBase.add( Integer.toString(i) );
		}
		
		
		
		int numSubstrings = 2000;
		HashSet<SequenceOfSymbols> substrings = HeuristicsForPickingBase.getBestSubstrings2(m, numSubstrings);
		System.out.println("Done choosing substrings");
		System.out.println(substrings);
		System.out.println();
		
		/*HashSet<SequenceOfSymbols> substrings = new HashSet<SequenceOfSymbols>();
		for (SequenceOfSymbols seq : m.keySet()) {
			substrings.addAll(seq.getSubstrings());
		}
		*/
		for (SequenceOfSymbols seq : m.keySet()) {
			currentBestDecomposition.put(seq, seq.getRawSequence().length());
		}
		
		System.out.println("Base size");
		System.out.println(maxBaseSize);
		System.out.println();
		System.out.println("Number of observations");
		System.out.println(m.size());
		System.out.println();
		System.out.println("Number of substrings: ");
		System.out.println(substrings.size());
		System.out.println();
		
		//int randomSample = 10;
		///HashSet<SequenceOfSymbols> sub = randomSubset(m, randomSample);
		
		while(currentBase.size() < maxBaseSize && substrings.size() > 0){
			PriorityQueue<SymbolCountPair> pq = new PriorityQueue<SymbolCountPair>();
			int i = 0;
			for (SequenceOfSymbols s : substrings) {
			
				HashSet<String> tempBase = (HashSet<String>) currentBase.clone();
				tempBase.add(s.getRawSequence());
				int improvement = 0;
								
				for (SequenceOfSymbols seq : m.keySet() ) {
					StringIntPair p;
					try {
						p = computeOptimalCompositionsNeccesaryCoinStyle(tempBase, seq.getRawSequence(), "Min");
					} catch (Exception e) {
						System.out.println(seq.getRawSequence());
						//System.out.println();
						e.printStackTrace();
						System.out.println("Problem computing best composition");
						return null;
					}
					 int extra = currentBestDecomposition.get(seq) - p.getI();
					 if (extra < 0){
						 System.out.println(extra);
						 System.out.println("Improvement less than 0 --> buggy");
						return null;
					 } 
					 improvement += extra*m.get(seq);
				} 
				//i++;
				//System.out.println(i);

				SymbolCountPair sc = new SymbolCountPair(-1*improvement, s);
				pq.add(sc);
			}
			SymbolCountPair bestAddition = pq.peek();
			substrings.remove(bestAddition.getSequence());
			currentBase.add(bestAddition.getSequence().getRawSequence());
			//System.out.println("Added to base:");
			//System.out.println(bestAddition.getSequence());
			//System.out.println();
			currentBestDecomposition = updateCurrentBestDecomposition(m.keySet(), currentBase);
		}
		
		HashSet<SequenceOfSymbols> returnBase = stringHashSetToSeqOfSymbols(currentBase);
		System.out.println(returnBase);
		return returnBase;
		
	}
	
	private static HashSet<SequenceOfSymbols> randomSubset(Map<SequenceOfSymbols, Integer> m, int randomSample) {
		HashSet<SequenceOfSymbols> sample = new HashSet<SequenceOfSymbols>();
		Random r = new Random();
		int count = 0;
		for (SequenceOfSymbols sequenceOfSymbols : m.keySet()) {
			count += m.get(sequenceOfSymbols);
		}
		
		for (int i = 0; i < randomSample; i++) {
		
			int rI = r.nextInt(count);
			int t = 0;
			for (SequenceOfSymbols sequenceOfSymbols : m.keySet()) {
				if (t>rI){
					sample.add(sequenceOfSymbols);
					break;
				}
				t+= m.get(sequenceOfSymbols);
				if (t==rI){
					sample.add(sequenceOfSymbols);
					break;
				}
			}
			
		}
		
		return sample;
	}
	private static HashMap<SequenceOfSymbols, Integer> updateCurrentBestDecomposition(Set<SequenceOfSymbols> set, HashSet<String> updatedBase){
		HashMap<SequenceOfSymbols, Integer> currentBestDecomp = new HashMap<SequenceOfSymbols, Integer>();
		for (SequenceOfSymbols sub : set) {
			StringIntPair opt;
			try {
				opt = computeOptimalCompositionsNeccesaryCoinStyle(updatedBase, sub.getRawSequence(), "Min");
				currentBestDecomp.put(sub, opt.getI());
			} catch (Exception e) {
				System.out.println("Problem computing optimal when updating currentBestDecomposition");
				e.printStackTrace();
				return null;
			}
		}
		return currentBestDecomp;
	}
	
	private static HashSet<SequenceOfSymbols> stringHashSetToSeqOfSymbols(HashSet<String> h){
		HashSet<SequenceOfSymbols> hr = new HashSet<SequenceOfSymbols>();
		for (String s : h) {
			hr.add(SequenceOfSymbols.fullStringToCompressed(s));
		}
		return hr;
	}
	
	public static HashSet<String> seqHashSetToStringOfSymbols(Set<SequenceOfSymbols> h){
		HashSet<String> hr = new HashSet<String>();
		for (SequenceOfSymbols s : h) {
			hr.add( s.getRawSequence() );
		}
		return hr;
	}
	
	private static SequenceOfSymbols getMinFromHashSet(Map<SequenceOfSymbols, Integer> scores) {
		double min = 0;
		SequenceOfSymbols minKey = null;
		boolean b = false;
		for (SequenceOfSymbols s : scores.keySet() ){
			if (b == false || scores.get(s) < min){
				min = scores.get(s);
				minKey = s;
				b = true;
			}
		}
		return minKey;
	}
	

	private static HashSet<SequenceOfSymbols> getBestSubstrings2(Map<SequenceOfSymbols, Integer> s, int numSubstrings){
		SymbolCounts sc = new SymbolCounts();
		for (SequenceOfSymbols sequenceOfSymbols : s.keySet()) {
			for (SequenceOfSymbols seqs : sequenceOfSymbols.getSubstrings()) {
				sc.updateFrequency(seqs, s.get(sequenceOfSymbols));
			}
		}
		
		PriorityQueue<SymbolCountPair> pq = new PriorityQueue<SymbolCountPair>();
		for (SequenceOfSymbols seq : sc.symbolToFrequency.keySet()) {
			int c = sc.symbolToFrequency.get(seq);
			SymbolCountPair scp = new SymbolCountPair(-1*c, seq);
			pq.add(scp);
		}
		
		HashSet<SequenceOfSymbols> r = new HashSet<SequenceOfSymbols>();
		for (int i = 0; i < numSubstrings; i++) {
			r.add( pq.remove().getSequence());
		}
		
		return r;
		
	}
	
	public static HashSet<SequenceOfSymbols> chooseTreeBase(int numDimensions, int length){
		HashSet<SequenceOfSymbols> h = new HashSet<SequenceOfSymbols>();

		if(length == 1){

			for (int i = 1; i <= numDimensions; i++) {
				h.add( SequenceOfSymbols.fullStringToCompressed(Integer.toString(i)) );
			}
		}
		else{

			for (int i = 1; i <= numDimensions; i++) {
				SequenceOfSymbols t = SequenceOfSymbols.fullStringToCompressed(Integer.toString(i) );
				HashSet<SequenceOfSymbols> rh = chooseTreeBase(numDimensions, length - 1);
				h.addAll(rh);
				for (SequenceOfSymbols sequenceOfSymbols : rh) {
					SequenceOfSymbols s = SequenceOfSymbols.concatenateSymbols(sequenceOfSymbols, t);
					h.add(s);
				}
				
			
			}
		}
		if (h.size() != ((int) (Math.pow(numDimensions, length+1) - 1 )/ (numDimensions - 1)) -1 ){
			System.out.println("Problem with tree method");
			System.out.println("Size of output doesn't seem right");
			System.out.println(((Math.pow(numDimensions, length+1) - 1 )/ (numDimensions - 1)) - 1 );
			System.out.println(h.size());
			System.out.println();

		}
		
		return h;

	}
	
	public static HashSet<SequenceOfSymbols> choosePowerBase(int numDimensions , int base, int maxExponent){
		HashSet<SequenceOfSymbols> h = new HashSet<SequenceOfSymbols>();

		for (int i = 1; i <= numDimensions; i++) {
			
			for (int j = 0; j <= maxExponent; j++) {
				SequenceOfSymbols s = new SequenceOfSymbols( i + ":" + Math.pow(base, j) );
				h.add(s);
			}
		}
		return h;

	}
	
	
	

}
