package hmm_sim;

import Jama.Matrix;

public class OutputDataPair {

	public Matrix getData() {
		return data;
	}
	public Matrix getxAxis() {
		return xAxis;
	}
	public OutputDataPair(Matrix data, Matrix xAxis) {
		super();
		this.data = data;
		this.xAxis = xAxis;
	}
	public Matrix data;
	public Matrix xAxis;
	
	
}
