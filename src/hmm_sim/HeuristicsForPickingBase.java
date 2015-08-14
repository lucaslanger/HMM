package hmm_sim;

import java.util.HashMap;
import java.util.HashSet;
import java.util.PriorityQueue;

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
		HashSet<SequenceOfSymbols> a = computeBaseSystem(seqs, 1, 4);
		System.out.println(a);
	}
	
	
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

	//Dynamic Programming algorithm
	public static int computeNumberCompositionsNeccesary(HashSet<String> base, String s){
		int[][] compositionsNeccesary = new int[s.length()][s.length()];
		try {
			return recursiveComputationsNeccessary(compositionsNeccesary, base, s, 0, s.length()-1);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Problem with computation of how to chop into substrings");
			return 0;
		}
		
	}
	
	private static int recursiveComputationsNeccessary(int[][] compositionsNeccesary, HashSet<String> base, String s, int i, int j) throws Exception{
		if(i==j){
			if(base.contains(  Character.toString(s.charAt(i)) ) ){
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
			//System.out.println("Memoization working");
			return compositionsNeccesary[i][j];
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
			compositionsNeccesary[i][j] = min;
			return min;
		}
	}
	
	public  static int formulaForStringScore(int compositionsNeccesary, int appearences){
		//MISSING THE LENGTH OF THE STRING
		return appearences - compositionsNeccesary;
	}

}
