package hmm_sim;

import java.util.HashMap;

public class SymbolInfo {
	
	private HashMap<String, Double> symbolToProbability;
	private HashMap<String, Integer> symbolToIndex;

	public SymbolInfo(HashMap<String, Double> symbolToProbability,
			HashMap<String, Integer> symbolToIndex) {
		super();
		this.symbolToProbability = symbolToProbability;
		this.symbolToIndex = symbolToIndex;
	}
	
	public HashMap<String, Double> getSymbolToProbability() {
		return symbolToProbability;
	}
	public HashMap<String, Integer> getSymbolToIndex() {
		return symbolToIndex;
	}
	


}
