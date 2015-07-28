package hmm_sim;

public class SymbolCountPair{
	
	public int getCount() {
		return count;
	}

	public String getSymbol() {
		return symbol;
	}

	private int count;
	private String symbol;

	public SymbolCountPair(int count, String symbol) {
		this.count = count;
		this.symbol = symbol;
	}

}
