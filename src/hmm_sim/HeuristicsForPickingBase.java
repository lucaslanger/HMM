package hmm_sim;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;

import javax.naming.BinaryRefAddr;

public class HeuristicsForPickingBase {
	
	public static void main(String[] args){
		String[] t = new String[]{"","","","",
				"1:7", "1:7", "1:7", "1:7",
				"1:14", "1:14", "1:14", "1:14",
				"1:21","1:21",
				"1:28",
				"1:35"
		};
		SequenceOfSymbols[] seqs = new SequenceOfSymbols[t.length];
		int i = 0;
		for(String s: t){
			seqs[i] = new SequenceOfSymbols(s);
			i++;
		}
		
		int numDimensions = 1;
		
		HashSet<String> base = new HashSet<String>();
		base.add("1");
		
		for (SequenceOfSymbols s: seqs) {
			computeOptimalCompositionsNeccesary(base, s.getRawSequence() );
		}
	}
	

	//Dynamic Programming algorithm
	public static int computeOptimalCompositionsNeccesary(HashSet<String> base, String s){
		HashMap<Integer, ArrayList<String>> possiblePlugins = new HashMap<Integer, ArrayList<String>>();
		for (int i = 0; i < s.length(); i++) {
			for (String string : base) {
				if (i + string.length() <  s.length() && string.equals(s.substring(i, i + string.length()) )){
					ArrayList<String> t = possiblePlugins.get(i);
					t.add(string);
					possiblePlugins.put(i, t);
				}
			}
		}
		
		int[][] compositionsNeccesary = new int[s.length()][s.length()];
		HashMap<String, String > bestCompositions = new HashMap<String, String>();
		try {
			return recursiveComputationsNeccessary(possiblePlugins, compositionsNeccesary, bestCompositions, base, s, 0, s.length()-1);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Problem with computation of how to chop into substrings");
			return 0;
		}
		
	}
	
	private static int recursiveComputationsNeccessary(HashMap<Integer, ArrayList<String>> possiblePlugins, int[][] compositionsNeccesary, HashMap<String, String> bestCompositions, HashSet<String> base, String s, int i, int j) throws Exception{
		if(i==j){
			if(base.contains(  Character.toString(s.charAt(i)) ) ){
				String s1 = Integer.toString(i) + ":" +  j;
				bestCompositions.put(s1, Character.toString(s1.charAt(i)) );
				return 1;
			}
			else{
				System.out.println("Basic Symbol not included!" );
				System.out.println( s.charAt(i) );
				System.out.println(base);
				throw new Exception();
			}
		}
		else if(compositionsNeccesary[i][j] != 0 ){

			return compositionsNeccesary[i][j];
		}
		else{
			int min = 0;
			int minIndex = 0;
			String minBaseString = "";
			boolean b = false;
			for (int h = i+1; h < j; h++) {
				for (String string : possiblePlugins.get(h)) {
					int cL = recursiveComputationsNeccessary(possiblePlugins, compositionsNeccesary, bestCompositions, base, s, i, h);
					int cR = recursiveComputationsNeccessary(possiblePlugins, compositionsNeccesary, bestCompositions, base, s, h+1+string.length(), j);
					int sum = cL + 1 +  cR;
					if (b == false || sum < min){
						min = sum;
						minIndex = h;
						minBaseString = string;
						b = true;
					}
				}
			}
			compositionsNeccesary[i][j] = min;
			
			String leftId = Integer.toString(i) +  ":" + minIndex;
			String rightId = Integer.toString(minIndex + minBaseString.length() ) + ":" + j;
			String lBest = bestCompositions.get(leftId);
			String rBest = bestCompositions.get(rightId);
			
			String combinedBest = lBest + minBaseString + rBest;
			
			String s1 = Integer.toString(i) + j;
			bestCompositions.put(s1, combinedBest);
			return min;
		}
	}
	
	public  static int formulaForStringScore(int compositionsNeccesary, int appearences){
		//MISSING THE LENGTH OF THE STRING
		return appearences - compositionsNeccesary;
	}
	
	
	/*
	public static SymbolCounts countSubstringOccurances(int numDimensions, SequenceOfSymbols[] strings){
		SymbolCounts substringOccurances = new SymbolCounts(numDimensions);
		for (SequenceOfSymbols sequenceOfSymbols : strings) {
			SymbolCounts t = sequenceOfSymbols.getSubstringCounts();
			for (SequenceOfSymbols substring : t.descKeySet() ) {
				substringOccurances.updateFrequency(substring, t.getSymbolToFrequency().get(substring) );
			}
		}
		return substringOccurances;
	}
	
	
	public static HashSet<SequenceOfSymbols> computeBaseSystemOld(SequenceOfSymbols[] strings, int numDimensions, int numberOfElementsInBase){
		HashSet<String > base = new HashSet<String>();
		for (int i = 1; i <= numDimensions; i++) {
			base.add( Integer.toString(i) );
			base.add( "" );
		}
		
		SymbolCounts subStringOccurances = countSubstringOccurances(numDimensions, strings);
		
		
		for (int i = 0; i < numberOfElementsInBase; i++) {
			HashMap<SequenceOfSymbols, Integer> scores = new HashMap<SequenceOfSymbols, Integer>();
			for (SequenceOfSymbols s : subStringOccurances.descKeySet()) {
				if (s.getSequence().equals("") == false && base.contains(s.getRawSequence()) == false){	// obviously dont want empty transition
					int cN = computeNumberCompositionsNeccesary(base, s.getRawSequence() );
					int occ = subStringOccurances.getSymbolToFrequency().get(s);
					scores.put(s, formulaForStringScore(cN, occ) );
				}
			}
			SequenceOfSymbols min = getMinFromHashSet(scores);
			base.add( min.getRawSequence() );
		}		
		
		HashSet<SequenceOfSymbols> rBase = new HashSet<SequenceOfSymbols>();
	
		for (String s : base) {
			SequenceOfSymbols t = SequenceOfSymbols.fullStringToCompressed(s);
			rBase.add(t);
		}
		
		return rBase;
		
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
	*/

}
