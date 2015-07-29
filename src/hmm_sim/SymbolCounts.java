package hmm_sim;

import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.TreeMap;

public class SymbolCounts {
	
	public TreeMap<String, Integer> getSymbolToFrequency() {
		return symbolToFrequency;
	}

	public int getDataCount() {
		return dataCount;
	}

	TreeMap<String, Integer> symbolToFrequency;
	int dataCount;
	
	public SymbolCounts(int numDimensions){
		this.symbolToFrequency = new TreeMap<String, Integer>();
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
	
	public NavigableMap<String, Integer> SortedKeys(){
		return symbolToFrequency.descendingMap();
	}

	public NavigableSet<String> descKeySet() {
		return symbolToFrequency.descendingKeySet();
	}

}
