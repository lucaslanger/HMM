package hmm_sim;

import java.util.HashMap;

public class SymbolCounts {
	
	HashMap<String, Integer> symbolToFrequency;
	int dataCount;
	public SymbolCounts(int numDimensions){
		this.symbolToFrequency = new HashMap<String, Integer>();
		this.dataCount = 0;
	}
	
	public void updateFrequency(String s){
		this.dataCount += 1;
		if (symbolToFrequency.containsKey(s)){
			symbolToFrequency.put(s, symbolToFrequency.get(s) + 1);
		}
		else{
			symbolToFrequency.put(s, 1);
		}
	}
	
	public SymbolInfo getProbabilityMap(){
		HashMap<String, Double> symbolToProbability = new HashMap<String, Double>();
		HashMap<String, Integer> symbolToIndex = new HashMap<String, Integer>();
		int c = 0;
		for (java.lang.String s : this.symbolToFrequency.keySet()) {
			symbolToProbability.put(s, ((double) this.symbolToFrequency.get(s))/ dataCount );
			symbolToIndex.put(s,  c);
		}
		return new SymbolInfo( symbolToProbability,  symbolToIndex );
	}

}
