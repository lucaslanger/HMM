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
	
	public ErrorPair searchOverBase(double[] mS, double[] bases, int maxK){
		String workingFolder = "keySearchPacMan/";
		String empModels = "Models_Emperical_" + workingFolder;
		String pltFolder = workingFolder + "Plotting_" + empModels + "/";
		
		FlowControl.createFolder(workingFolder);

		LabyrinthGraph l = LabyrinthGraph.pacMan(workingFolder, hSize, stretchFactor, key, false);
		l.generateData(trajectorySizes, repetitions, false);
		
		FlowControl.readDataIntoModels(workingFolder, basisSize);

		double[][] modelSizes = new double[][]{mS};
		double[][] basesToTest = new double[][]{bases};
		double[][] errorTestingVSModelSize = new double[basesToTest[0].length][modelSizes[0].length];
		double[][] errorTrainingVSModelSize = new double[basesToTest[0].length][modelSizes[0].length];

		double[][] xaxis = new double[basesToTest[0].length][modelSizes[0].length];

		int[] shortestPaths = l.shortestPathsFromKey();

		System.out.println("HACK ON PSEUDOINVERSE");
		System.out.println("HACK ON PSEUDOINVERSE");
		System.out.println("HACK ON PSEUDOINVERSE");
			
		ModelRetrieval mr = new ModelRetrieval(workingFolder, empModels, "Models_True_" + workingFolder, basisSize, base);
		ModelRetrieval trueModel = new ModelRetrieval(workingFolder, empModels, "Models_True_" + workingFolder, basisSize, base);
		HankelSVDModel trueHSVD = trueModel.readTrueModel("Models_True_" + workingFolder);
		int rank = trueHSVD.getRank();
		System.out.println("Rank");
		System.out.println(rank);
		QueryEngine trueQE = trueHSVD.buildHankelBasedModel(base, rank); 
		
		QueryEngine[] learnedModels = new QueryEngine[modelSizes[0].length];
		for (int i = 0; i < learnedModels.length; i++) {
			learnedModels[i] = mr.getSpecificModelSizeQueryEngines(1, (int) modelSizes[0][i]).get(dataSizeForFixedPlots)[0];
			checkQEForDifferences(trueQE, learnedModels[i], bases, 500);
		}
		
		Matrix[][][] alphaKStatesOver_MS_bases = new Matrix[modelSizes[0].length][basesToTest[0].length][maxK];
		for (int i = 0; i < alphaKStatesOver_MS_bases.length; i++) {
			for (int j = 0; j < basesToTest[0].length; j++) {
				alphaKStatesOver_MS_bases[i][j] = learnedModels[i].getAllKStateQueries(maxK, (int) basesToTest[0][j], base);
			}
		}

		for (int r = 0; r < repetitions; r++) {
			
			HashMap<String, int[]> trainingSamples = l.createObservationDistanceSamples(shortestPaths, maxK, samples);
			HashMap<String, int[]> testingSamples = l.createObservationDistanceSamples(shortestPaths, maxK, samples);
			
			for (int i = 0; i < modelSizes[0].length; i++) {
				
				//QueryEngine learnedModel = learnedModels[i];
				
				for (int j = 0; j < basesToTest[0].length; j++) {
				
					int maxPow = (int) basesToTest[0][j];
					//Matrix[] alphaKStates = learnedModel.getAllKStateQueries(maxK, maxPow, base);
					Matrix[] alphaKStates = alphaKStatesOver_MS_bases[i][j];
					
					Matrix theta = l.getAlphaFromSampledData(trainingSamples, alphaKStates);
		
					//System.out.print("Base: ");
					//System.out.println(maxPow);
					double eTraining = KeySearching.determineError(theta, alphaKStates, trainingSamples);
					//System.out.println(eTraining);
					
					double eTesting = KeySearching.determineError(theta, alphaKStates, testingSamples );
					//System.out.println(eTesting);
					
					errorTestingVSModelSize[j][i] += eTesting;
					errorTrainingVSModelSize[j][i] += eTraining;
				}
				
				
			}
		}
		for (int j = 0; j < basesToTest[0].length; j++) {
			for (int i = 0; i < modelSizes[0].length; i++) {
				xaxis[j][i] = modelSizes[0][i];
				errorTestingVSModelSize[j][i] /= repetitions*samples;
				errorTrainingVSModelSize[j][i] /= repetitions*samples;
			}
		}
		
		Matrix errTraining = new Matrix(errorTrainingVSModelSize);
		Matrix errTesting = new Matrix(errorTestingVSModelSize);
		
		errTraining.print(5, 5);
		errTesting.print(5, 5);
		
		String title = "Distance Predictions";
		String internalComment = "Darker Curves --> Richer Base System";
		OutputData.outputData(pltFolder + "KeyFindingErrorTraining_MaxK:" + maxK , "Number Of States", "Error Norm_2", xaxis, errorTrainingVSModelSize, title, internalComment);
		OutputData.outputData(pltFolder + "KeyFindingErrorTesting_MaxK:" + maxK, "Number Of States", "Error Norm_2", xaxis, errorTestingVSModelSize, title, internalComment);

		return new ErrorPair(errorTrainingVSModelSize, errorTestingVSModelSize);
	}
	
	private void checkQEForDifferences(QueryEngine trueModel, QueryEngine learnedModel, double[] baseLimits, int maxquery) {
		double[][] probabilities = new double[baseLimits.length][maxquery];
		double[] differences = new double[baseLimits.length];
		
		double[] tProb = new double[maxquery];
		for (int i = 0; i < tProb.length; i++) {
			tProb[i] = trueModel.probabilityQuery(i, 64, base, true);
		}
		
		for (int i = 0; i < probabilities.length; i++) {
			int m = (int) baseLimits[i];
			for (int j = 0; j < maxquery; j++) {
				probabilities[i][j] = learnedModel.probabilityQuery(j, m, base, true);
			}
		}
		
		for (int i = 0; i < probabilities.length; i++) {
			for (int j = 0; j < maxquery; j++) {
				differences[i] += Math.abs(tProb[j] - probabilities[i][j]);
			}
		}
		
		System.out.println("Differences");
		System.out.println( Arrays.toString(differences) );
		
	}

	public static double determineError(Matrix theta, Matrix[] alphaKStates, HashMap<String, int[]> testSamples){
		int widthOfA = alphaKStates[0].getArrayCopy().length;
		int numSamples = testSamples.get("Durations").length;
		double[][] a = new double[numSamples][widthOfA];
		double[][] b = new double[numSamples][1];
 		for (int i = 0; i < numSamples; i++) {
			int traj = testSamples.get("Durations")[i];
			a[i] = alphaKStates[traj].getArrayCopy()[0];
			b[i][0] = testSamples.get("Distances")[i];
		}
		
 		Matrix A = new Matrix(a);
 		
 		Matrix B = new Matrix(b);
 		
 		/*System.out.println(A.times(theta).minus(B).norm2());
 		System.out.println( norm2Custom(A.times(theta).minus(B).transpose().getArray()[0]) );
 		System.out.println(A.times(theta).minus(B).norm1() );
 		System.out.println( norm1Custom(A.times(theta).minus(B).transpose().getArray()[0]) );

 		System.out.println();
 		*/
		return A.times(theta).minus(B).norm2();
	}

	public ErrorPair searchOverMaxK(double[] mS, double[] maxKs, int maxPow){
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

		System.out.println("HACK ON PSEUDOINVERSE");
		System.out.println("HACK ON PSEUDOINVERSE");
		System.out.println("HACK ON PSEUDOINVERSE");
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
					/*int trueModelRank = mr.readTrueModel("Models_True_" + workingFolder).getRank();
					System.out.println("TrueModel rank: ");
					System.out.println(trueModelRank);
					System.out.println("Rep: " + r + " MS: " + m + " MaxK: " + k);
					System.out.println();
					*/
					QueryEngine learnedModel = mr.getSpecificModelSizeQueryEngines(repetitions, m).get(dataSizeForFixedPlots)[0];
					Matrix[] alphaKStates = learnedModel.getAllKStateQueries(k, maxPow, base);
					
					Matrix theta = l.getAlphaFromSampledData(trainingSamples, alphaKStates);
		
					double eTesting = KeySearching.determineError(theta, alphaKStates, testingSamples  );
					double eTraining = KeySearching.determineError(theta, alphaKStates, trainingSamples );
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
		
		String title = "Distance Predictions";
		String internalComment = "Darker Curves --> Richer Base System";
		OutputData.outputData(pltFolder + "KeyFindingErrorTraining_Base:" + base , "Number Of Statesze", "Error Norm2()", xaxis, errorTrainingVSModelSize, title, internalComment);
		OutputData.outputData(pltFolder + "KeyFindingErrorTesting_Base:" + base, "Number Of States", "Error Norm2()", xaxis, errorTestingVSModelSize, title, internalComment);

		return new ErrorPair(errorTrainingVSModelSize, errorTestingVSModelSize);
	}

}
