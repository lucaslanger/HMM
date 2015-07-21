package hmm_sim;

public class IntegerPair {
	
	private int lowerInt;
	private int upperInt;

	
	public IntegerPair(int i1, int i2) {
		if (i1 > i2){
			this.upperInt = i1;
			this.lowerInt = i2;
		}
		else{
			this.upperInt = i1;
			this.lowerInt = i2;
		}
	}

	public int getLowerInt() {
		return lowerInt;
	}

	public int getUpperInt() {
		return upperInt;
	}

}
