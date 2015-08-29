package hmm_sim;

import Jama.Matrix;

public class OutputDataPair {

	public Matrix getData() {
		return data;
	}
	public Matrix getxAxis() {
		return xAxis;
	}
	public Matrix getSpreads(){
		return SPREADS;
	}
	
	public OutputDataPair(Matrix data, Matrix xAxis, Matrix SPREADS) {
		super();
		this.data = data;
		this.xAxis = xAxis;
		this.SPREADS = SPREADS;
	}
	public Matrix data;
	public Matrix xAxis;
	public Matrix SPREADS;
	
	
}
