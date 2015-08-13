package hmm_sim;

import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

public class HeuristicsForPickingBase {
	
	
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
	
	public static HashSet<SequenceOfSymbols> computeBaseSystem(SequenceOfSymbols[] strings, int numDimensions, int numberOfElementsInBase){
		HashSet<String > base = new HashSet<String>();
		for (int i = 1; i <= numDimensions; i++) {
			base.add( Integer.toString(i) );
		}
		
		SymbolCounts subStringOccurances = countSubstringOccurances(numDimensions, strings);
		
		HashMap<SequenceOfSymbols, Integer> scores = new HashMap<SequenceOfSymbols, Integer>();
		
		for (int i = 0; i < numberOfElementsInBase; i++) {
			for (SequenceOfSymbols s : subStringOccurances.descKeySet()) {
				int cN = computeNumberCompositionsNeccesary(base, s.getRawSequence() );
				int occ = subStringOccurances.getSymbolToFrequency().get(s);
				scores.put(s, formulaForStringScore(cN, occ) );
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

	//Dynamic Programming algorithm
	public static int computeNumberCompositionsNeccesary(HashSet<String> base, String s){
		int[][] compositionsNeccesary = new int[s.length()][s.length()];
		try {
			return recursiveComputationsNeccessary(compositionsNeccesary, base, s, 0, s.length()-1);
		} catch (Exception e) {
			System.out.println("Problem with computation of how to chop into substrings");
			return 0;
		}
		
	}
	
	private static int recursiveComputationsNeccessary(int[][] compositionsNeccesary, HashSet<String> base, String s, int i, int j) throws Exception{
		if(i==j){
			if(base.contains(s.charAt(i)) ){
				return 1;
			}
			else{
				System.out.println("Basic Symbol not included!" );
				throw new Exception();
			}
		}
		else{
			int min = 0;
			boolean b = false;
			for (int h = i+1; h < j; h++) {
				int cL = recursiveComputationsNeccessary(compositionsNeccesary, base, s, i, h);
				int cR = recursiveComputationsNeccessary(compositionsNeccesary, base, s, h+1, j);
				int sum = cL + cR;
				if (b == false || sum < min){
					min = sum;
					b = true;
				}
			}
			return min;
		}
	}
	
	public  static int formulaForStringScore(int compositionsNeccesary, int appearences){
		return appearences - compositionsNeccesary;
	}

}
