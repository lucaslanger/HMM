package hmm_sim;

import java.util.Arrays;
import java.util.HashMap;

import Jama.Matrix;

public class KeySearching {
	
	
	private int samples;
	private int key;
	private int basisSize;
	private int hSize;
	private int stretchFactor;
	private int repetitions;
	private int base;
	private int[] trajectorySizes;
	private int dataSizeForFixedPlots;
	
	public KeySearching(int samples, int key, int basisSize, int hSize, int stretchFactor, int[] trajectorySizes, int dataSizeForFixedPlots, int repetitions, int base) {
		this.samples = samples;
		this.key = key;
		this.basisSize = basisSize;
		this.hSize = hSize;
		this.stretchFactor = stretchFactor;
		this.repetitions = repetitions;
		this.base = base;
		this.trajectorySizes = trajectorySizes;
		this.dataSizeForFixedPlots = dataSizeForFixedPlots;
	}

	public void search(double[] mS, double[] maxKs){
		String workingFolder = "keySearchPacMan/";
		String empModels = "Models_Emperical_" + workingFolder;
		String pltFolder = workingFolder + "Plotting_" + empModels + "/";
		
		FlowControl.createFolder(workingFolder);

		LabyrinthGraph l = LabyrinthGraph.pacMan(workingFolder, hSize, stretchFactor, key);
		l.generateData(trajectorySizes, repetitions, false);
		
		FlowControl.readDataIntoModels(workingFolder, basisSize);

		double[][] modelSizes = new double[][]{mS};
		double[][] maxKsToTest = new double[][]{maxKs};
		double[][] errorTestingVSModelSize = new double[maxKsToTest[0].length][modelSizes[0].length];
		double[][] errorTrainingVSModelSize = new double[maxKsToTest[0].length][modelSizes[0].length];

		double[][] xaxis = new double[maxKsToTest[0].length][modelSizes[0].length];

		int[] shortestPaths = l.shortestPathsFromKey();

		System.out.println("AVERAGING");
		System.out.println("HACK ON PSEUDOINVERSE");
			
		for (int j = 0; j < maxKsToTest[0].length; j++) {
			int k = (int) maxKsToTest[0][j]; 
			//double[][] trueDistanceAhead = l.dynamicallyDetermineTrueDistanceKAhead(shortestPaths, k);
			for (int r = 0; r < repetitions; r++) {
				HashMap<String, int[]> trainingSamples = l.createObservationDistanceSamples(shortestPaths, k, samples);
				HashMap<String, int[]> testingSamples = l.createObservationDistanceSamples(shortestPaths, k, samples);

				for (int i = 0; i < modelSizes[0].length; i++) {
					int m = (int) modelSizes[0][i];
					ModelRetrieval mr = new ModelRetrieval(workingFolder, empModels, "Models_True_" + workingFolder, basisSize, base);
					
					System.out.println("Rep: " + r + " MS: " + m + " MaxK: " + k);
				
					QueryEngine learnedModel = mr.getSpecificModelSizeQueryEngines(repetitions, m).get(dataSizeForFixedPlots)[0];
					Matrix[] alphaKStates = learnedModel.getAllKStateQueries(k, base);
					
					Matrix theta = l.getAlphaFromSampledData(trainingSamples, alphaKStates);
		
					double eTesting = l.determineError(theta, alphaKStates, testingSamples);
					double eTraining = l.determineError(theta, alphaKStates, trainingSamples);
					errorTestingVSModelSize[j][i] += eTesting;
					errorTrainingVSModelSize[j][i] += eTraining;
					
				}
				
			}
			
		}
		for (int j = 0; j < maxKsToTest[0].length; j++) {
			for (int i = 0; i < modelSizes[0].length; i++) {
				xaxis[j][i] = modelSizes[0][i];
				errorTestingVSModelSize[j][i] /= repetitions;
				errorTrainingVSModelSize[j][i] /= repetitions;
			}
		}
	
		
		Matrix errTraining = new Matrix(errorTrainingVSModelSize);
		Matrix errTesting = new Matrix(errorTestingVSModelSize);
		System.out.println("Training error:");
		errTraining.print(5, 5);
		System.out.println("Testing error");
		errTesting.print(5, 5);
		
		OutputData.outputData(pltFolder + "KeyFindingErrorTraining", "ModelSize | NOTE: Lighter curves --> Lower Trajectory Lengths", "Error Norm1()", xaxis, errorTrainingVSModelSize);
		OutputData.outputData(pltFolder + "KeyFindingErrorTesting", "ModelSize | NOTE: Lighter curves --> Lower Trajectory Lengths", "Error Norm1()", xaxis, errorTestingVSModelSize);

	}

}
