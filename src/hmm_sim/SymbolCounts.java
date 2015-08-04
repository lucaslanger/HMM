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
	
	public int getNumDimensions(){
		return this.numDimensions;
	}

	TreeMap<SequenceOfSymbols, Integer> symbolToFrequency;
	int dataCount;
	private int numDimensions;
	
	public SymbolCounts(int numDimensions){
		this.symbolToFrequency = new TreeMap<SequenceOfSymbols, Integer>();
		this.dataCount = 0;
		this.numDimensions = numDimensions;
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
	
	public void insertKeyTreatAsHashSet(SequenceOfSymbols s){
		if (symbolToFrequency.containsKey(s) == false){
			symbolToFrequency.put(s, 1);
			this.dataCount += 1 ;
		}
	}

	public NavigableSet<SequenceOfSymbols> descKeySet() {
		//System.out.println("Double check that this sorts by value not key");
		return symbolToFrequency.descendingKeySet();
	}
	
	public NavigableSet<SequenceOfSymbols> incrKeySet() {
		//System.out.println("Double check that this sorts by value not key");
		return symbolToFrequency.navigableKeySet();
	}
	
	public String toString(){
		String output = "";
		for (SequenceOfSymbols s : this.symbolToFrequency.keySet()) {
			output = output + "Sequence: " + s.getSequence() + " Occurences: " + this.symbolToFrequency.get(s) + ", "; 
		}
		return output;
	}

}
