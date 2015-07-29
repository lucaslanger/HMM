package hmm_sim;

import java.util.HashMap;

public class SymbolCounts {
	
	public HashMap<String, Integer> getSymbolToFrequency() {
		return symbolToFrequency;
	}

	public int getDataCount() {
		return dataCount;
	}

	HashMap<String, Integer> symbolToFrequency;
	int dataCount;
	public SymbolCounts(int numDimensions){
		this.symbolToFrequency = new HashMap<String, Integer>();
		this.dataCount = 0;
	}
	
	public void updateFrequency(String s, int update){
		this.dataCount += update;
		if (symbolToFrequency.containsKey(s)){
			symbolToFrequency.put(s, symbolToFrequency.get(s) + update);
		}
		else{
			symbolToFrequency.put(s, update);
		}
	}

}
