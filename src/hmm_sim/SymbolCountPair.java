package hmm_sim;

public class SymbolCountPair implements Comparable{
	
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

	@Override
	public int compareTo(Object arg0) {
		SymbolCountPair p2 = (SymbolCountPair) arg0;
		if ( this.getCount() < p2.getCount() ){
			return -1;
		}
		else if( this.getCount() == p2.getCount()){
			return 0;
		}
		else{
			return 1;
		}
		
	}
}


