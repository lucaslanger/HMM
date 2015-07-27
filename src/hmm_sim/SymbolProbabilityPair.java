package hmm_sim;

public class SymbolProbabilityPair{
	
	public double getProbability() {
		return probability;
	}

	public String getSymbol() {
		return symbol;
	}

	private double probability;
	private String symbol;

	public SymbolProbabilityPair(double probability, String symbol) {
		super();
		this.probability = probability;
		this.symbol = symbol;
	}

}
