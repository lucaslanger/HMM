package hmm_sim;

public class SymbolCountPair{
	
	public int getCount() {
		return count;
	}

	public SequenceOfSymbols getSequence() {
		return sequence;
	}

	private int count;
	private SequenceOfSymbols sequence;

	public SymbolCountPair(int count, SequenceOfSymbols sequence) {
		this.count = count;
		this.sequence = sequence;
	}

}
