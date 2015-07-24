package hmm_sim;

public class ErrorPair {

	double[][] trainingErrors;
	double[][] testingErrors;

	public ErrorPair(double[][] trainingErrors, double[][] testingErrors) {
		super();
		this.trainingErrors = trainingErrors;
		this.testingErrors = testingErrors;
	}
	
	public double[][] getTrainingErrors() {
		return trainingErrors;
	}

	public double[][] getTestingErrors() {
		return testingErrors;
	}

}
