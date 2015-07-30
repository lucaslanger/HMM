package hmm_sim;

public class SymbolCountPair{
	
	public int getCount() {
		return count;
	}

	public SequenceOfSymbols getSymbol() {
		return symbol;
	}

	private int count;
	private SequenceOfSymbols symbol;

	public SymbolCountPair(int count, SequenceOfSymbols symbol) {
		this.count = count;
		this.symbol = symbol;
	}

}
