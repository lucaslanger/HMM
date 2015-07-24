package hmm_sim;

import java.util.Arrays;

import Jama.Matrix;

public class ModelKeySearchComparison {
	
	double[][][] trainingData;
	double[][][] testingData;
	double[][][] xAxes;
	public ModelKeySearchComparison(double[][][] trainingData,
			double[][][] testingData, double[][][] xAxes) {
		super();
		this.trainingData = trainingData;
		this.testingData = testingData;
		this.xAxes = xAxes;
	}
	public double[][][] getTrainingData() {
		return trainingData;
	}
	public double[][][] getTestingData() {
		return testingData;
	}
	public double[][][] getxAxes() {
		return xAxes;
	}
	
	public void printOut() {
		System.out.println("BlocksDown --> Increasing MaxK");
		System.out.println("Individual Block Format: --> Higher ModelSize , Downwards for Higher Base");
		System.out.println();
		for (int i = 0; i < testingData.length; i++) {
			Matrix train = new Matrix(trainingData[i]);
			Matrix test = new Matrix(testingData[i]);
			System.out.println("train");
			train.print(5, 5);
			System.out.println("test");
			test.print(5, 5);
		}
		System.out.println("Done");
		System.out.println();
	}
	
	

}
