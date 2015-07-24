package hmm_sim;

import java.io.FileInputStream;
import java.io.ObjectInputStream;
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
	private int maxPower;
	
	public KeySearching(int samples, int key, int basisSize, int hSize, int stretchFactor, int[] trajectorySizes, int dataSizeForFixedPlots, int repetitions, int maxPower, int base) {
		this.samples = samples;
		this.key = key;
		this.basisSize = basisSize;
		this.hSize = hSize;
		this.stretchFactor = stretchFactor;
		this.repetitions = repetitions;
		this.base = base;
		this.trajectorySizes = trajectorySizes;
		this.dataSizeForFixedPlots = dataSizeForFixedPlots;
		this.maxPower = maxPower;
	}

	public ErrorPair search(double[] mS, double[] maxKs){
		String workingFolder = "keySearchPacMan/";
		String empModels = "Models_Emperical_" + workingFolder;
		String pltFolder = workingFolder + "Plotting_" + empModels + "/";
		
		FlowControl.createFolder(workingFolder);

		LabyrinthGraph l = LabyrinthGraph.pacMan(workingFolder, hSize, stretchFactor, key, false);
		l.generateData(trajectorySizes, repetitions, false);
		
		FlowControl.readDataIntoModels(workingFolder, basisSize);

		double[][] modelSizes = new double[][]{mS};
		double[][] maxKsToTest = new double[][]{maxKs};
		double[][] errorTestingVSModelSize = new double[maxKsToTest[0].length][modelSizes[0].length];
		double[][] errorTrainingVSModelSize = new double[maxKsToTest[0].length][modelSizes[0].length];

		double[][] xaxis = new double[maxKsToTest[0].length][modelSizes[0].length];

		int[] shortestPaths = l.shortestPathsFromKey();

		//System.out.println("HACK ON PSEUDOINVERSE");
			
		for (int j = 0; j < maxKsToTest[0].length; j++) {
			int k = (int) maxKsToTest[0][j]; 
			//double[][] trueDistanceAhead = l.dynamicallyDetermineTrueDistanceKAhead(shortestPaths, k);
			for (int r = 0; r < repetitions; r++) {
				HashMap<String, int[]> trainingSamples = l.createObservationDistanceSamples(shortestPaths, k, samples);
				HashMap<String, int[]> testingSamples = l.createObservationDistanceSamples(shortestPaths, k, samples);

				for (int i = 0; i < modelSizes[0].length; i++) {
					int m = (int) modelSizes[0][i];
					ModelRetrieval mr = new ModelRetrieval(workingFolder, empModels, "Models_True_" + workingFolder, basisSize, base);
					/*int trueModelRank = mr.readTrueModel("Models_True_" + workingFolder).getRank();
					System.out.println("TrueModel rank: ");
					System.out.println(trueModelRank);
					System.out.println("Rep: " + r + " MS: " + m + " MaxK: " + k);
					System.out.println();
					*/
					QueryEngine learnedModel = mr.getSpecificModelSizeQueryEngines(repetitions, m).get(dataSizeForFixedPlots)[0];
					Matrix[] alphaKStates = learnedModel.getAllKStateQueries(k, this.maxPower, base);
					
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
				errorTestingVSModelSize[j][i] /= repetitions*samples;
				errorTrainingVSModelSize[j][i] /= repetitions*samples;
			}
		}
	
		
		Matrix errTraining = new Matrix(errorTrainingVSModelSize);
		Matrix errTesting = new Matrix(errorTestingVSModelSize);
		
		/*System.out.println("Training error:");
		errTraining.print(5, 5);
		System.out.println("Testing error");
		errTesting.print(5, 5);
		*/
		
		OutputData.outputData(pltFolder + "KeyFindingErrorTraining_Base:" + this.maxPower , "ModelSize | NOTE: Lighter curves --> Lower Trajectory Lengths", "Error Norm2()", xaxis, errorTrainingVSModelSize);
		OutputData.outputData(pltFolder + "KeyFindingErrorTesting_Base:" + this.maxPower, "ModelSize | NOTE: Lighter curves --> Lower Trajectory Lengths", "Error Norm2()", xaxis, errorTestingVSModelSize);

		return new ErrorPair(errorTrainingVSModelSize, errorTestingVSModelSize);
	}

}
