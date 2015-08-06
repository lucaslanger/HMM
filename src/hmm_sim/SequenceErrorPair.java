package hmm_sim;

public class SequenceErrorPair implements Comparable<SequenceErrorPair>{
	
	private SequenceOfSymbols seq;
	public SequenceErrorPair(SequenceOfSymbols seq, double error) {
		super();
		this.seq = seq;
		this.error = error;
	}
	public SequenceOfSymbols getSeq() {
		return seq;
	}
	public double getError() {
		return error;
	}
	private double error;
	
	@Override
	public int compareTo(SequenceErrorPair o) {
		if (this.getError() > o.getError()){
			return 1;
		}
		else if (this.getError() == o.getError()){
			return 0;
		}
		else{
			return -1;
		}
	}
	
	public String toString(){
		return "Sequence: " + this.getSeq() + ", " + "Error: " + -1*this.getError();
	}
	

	

}
