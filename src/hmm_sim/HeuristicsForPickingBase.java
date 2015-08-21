package hmm_sim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.PriorityQueue;
import java.util.Set;

import javax.naming.BinaryRefAddr;

public class HeuristicsForPickingBase {
	
	public static void main(String[] args){
		//System.out.println("dsfsf".substring(1, 3));
		HeuristicsForPickingBase.selectBaseTest();
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
		base.add("2212");
		//base.add("221");
		base.add("21");
				
		for (SequenceOfSymbols s: seqs) {
			try {
				StringIntPair opt = computeOptimalCompositionsNeccesary(base, s.getRawSequence(), "Max" );
				System.out.println(s);
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
	

	//Dynamic Programming algorithm
	public static StringIntPair computeOptimalCompositionsNeccesary(HashSet<String> base, String s, String type) throws Exception{
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
		
		HashSet<SequenceOfSymbols> substrings = new HashSet<SequenceOfSymbols>();
		for (SequenceOfSymbols seq : m.keySet()) {
			substrings.addAll(seq.getSubstrings());
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
		
		while(currentBase.size() < maxBaseSize && substrings.size() > 0){
			PriorityQueue<SymbolCountPair> pq = new PriorityQueue<SymbolCountPair>();
			int i = 0;
			for (SequenceOfSymbols s : substrings) {
				//System.out.println(s);
				//System.out.println(i);
				//i++;
				HashSet<String> tempBase = (HashSet<String>) currentBase.clone();
				tempBase.add(s.getRawSequence());
				int improvement = 0;
				for (SequenceOfSymbols seq : m.keySet()) {
					//System.out.println(seq);
					StringIntPair p;
					try {
						p = computeOptimalCompositionsNeccesary(tempBase, seq.getRawSequence(), "Min");
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

				SymbolCountPair sc = new SymbolCountPair(-1*improvement, s);
				pq.add(sc);
			}
			SymbolCountPair bestAddition = pq.peek();
			substrings.remove(bestAddition.getSequence());
			currentBase.add(bestAddition.getSequence().getRawSequence());
			//System.out.println("Added to base:");
			//System.out.println(bestAddition.getSequence());
			currentBestDecomposition = updateCurrentBestDecomposition(m.keySet(), currentBase);
		}
		
		HashSet<SequenceOfSymbols> returnBase = stringHashSetToSeqOfSymbols(currentBase);
		System.out.println(returnBase);
		return returnBase;
		
	}
	
	private static HashMap<SequenceOfSymbols, Integer> updateCurrentBestDecomposition(Set<SequenceOfSymbols> set, HashSet<String> updatedBase){
		HashMap<SequenceOfSymbols, Integer> currentBestDecomp = new HashMap<SequenceOfSymbols, Integer>();
		for (SequenceOfSymbols sub : set) {
			StringIntPair opt;
			try {
				opt = computeOptimalCompositionsNeccesary(updatedBase, sub.getRawSequence(), "Min");
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
	
	private static SequenceOfSymbols getMinFromHashSet(HashMap<SequenceOfSymbols, Integer> scores) {
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
	
	/*
	private static SequenceOfSymbols[] getBestSubstrings(SequenceOfSymbols[] s, int numSubstrings){
		SymbolCounts sc = new SymbolCounts();
		for (SequenceOfSymbols sequenceOfSymbols : s) {
		
		}
		
	}
	*/

}
