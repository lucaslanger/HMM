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
	
	public HashMap<String, Double>() getProbabilityMap(){
		HashMap<String, Double> symbolToProbability = new HashMap<String, Double>();
		for (String s : this.symbolToFrequency.keySet()) {
			symbolToProbability.put(s, ((double) this.symbolToFrequency.get(s))/ dataCount );
		}
		return symbolToProbability;
	}

}
