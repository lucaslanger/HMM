package hmm_sim;

import java.util.NavigableMap;
import java.util.NavigableSet;
import java.util.TreeMap;

public class SymbolCounts {
	
	public TreeMap<SequenceOfSymbols, Integer> getSymbolToFrequency() {
		return symbolToFrequency;
	}

	public int getDataCount() {
		return dataCount;
	}

	TreeMap<SequenceOfSymbols, Integer> symbolToFrequency;
	int dataCount;
	
	public SymbolCounts(int numDimensions){
		this.symbolToFrequency = new TreeMap<SequenceOfSymbols, Integer>();
		this.dataCount = 0;
	}
	
	public void updateFrequency(SequenceOfSymbols s, int update){
		this.dataCount += update;
		if (symbolToFrequency.containsKey(s)){
			symbolToFrequency.put(s, symbolToFrequency.get(s) + update);
		}
		else{
			symbolToFrequency.put(s, update);
		}
	}
	
	public NavigableMap<SequenceOfSymbols, Integer> SortedKeys(){
		return symbolToFrequency.descendingMap();
	}

	public NavigableSet<SequenceOfSymbols> descKeySet() {
		return symbolToFrequency.descendingKeySet();
	}
	
	public String toString(){
		String output = "";
		for (SequenceOfSymbols s : this.symbolToFrequency.keySet()) {
			output = output + "Sequence: " + s.getSequence() + " Occurences: " + this.symbolToFrequency.get(s) + ", "; 
		}
		return output;
	}

}
