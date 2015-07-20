package hmm_sim;

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

	public void search(){
		String workingFolder = "keySearchPacMan/";
		String empModels = "Models_Emperical_" + workingFolder;
		String pltFolder = workingFolder + "Plotting_" + empModels + "/";
		
		FlowControl.createFolder(workingFolder);

		LabyrinthGraph l = LabyrinthGraph.pacMan(workingFolder, hSize, stretchFactor, key);
		l.generateData(trajectorySizes, repetitions, false);
		
		FlowControl.readDataIntoModels(workingFolder, basisSize);

		double[][] modelSizes = new double[][]{{5, 10, 20 ,40, 80}};
		double[][] maxKsToTest = new double[][]{{5, 20, 40, 80, 160, 320, 640, 1280}};
		double[][] errorVSModelSize = new double[maxKsToTest[0].length][modelSizes[0].length];
		double[][] xaxis = new double[maxKsToTest[0].length][modelSizes[0].length];

		int[] shortestPaths = l.shortestPathsFromKey();

		System.out.println("AVERAGING");
		System.out.println("HACK ON PSEUDOINVERSE");
			
		for (int j = 0; j < maxKsToTest[0].length; j++) {
			int k = (int) maxKsToTest[0][j]; 
			HashMap<String, int[]> trainingSamples = l.createObservationDistanceSamples(shortestPaths, k, samples);
			HashMap<String, int[]> testingSamples = l.createObservationDistanceSamples(shortestPaths, k, samples);

			//double[][] trueDistanceAhead = l.dynamicallyDetermineTrueDistanceKAhead(shortestPaths, k);
			
			for (int i = 0; i < modelSizes[0].length; i++) {	
				for (int r = 0; r < repetitions; r++) {
						
					int m = (int) modelSizes[0][i];
					//testEngine a = new testEngine(workingFolder, empModels, "Models_True_" + workingFolder, dataSizeForFixedPlots , basisSize, base, new int[]{}, m, 1 , false);
					//QueryEngine learnedModel = a.fixedModelQE.get(dataSizeForFixedPlots)[0];
	
					ModelRetrieval mr = new ModelRetrieval(workingFolder, empModels, "Models_True_" + workingFolder, basisSize, base);
					QueryEngine learnedModel = mr.getSpecificModelSizeQueryEngines(repetitions, m).get(dataSizeForFixedPlots)[r];
					
					Matrix[] alphaKStates = learnedModel.getAllKStateQueries(k, base);
					
					Matrix theta = l.getAlphaFromSampledData(trainingSamples, alphaKStates);
		
					double e = l.determineError(theta, alphaKStates, testingSamples);
					errorVSModelSize[j][i] += e;
				}
				xaxis[j][i] = modelSizes[0][i];
				errorVSModelSize[j][i] /= repetitions;
			}
			
		}
		Matrix pe = new Matrix(errorVSModelSize);
		pe.print(5, 5);
		
		System.out.println("Writing data to: " + pltFolder + "KeyFindingError");
		OutputData.outputData(pltFolder + "KeyFindingError", "ModelSize | NOTE: Lighter curves --> Lower MaxK", "Error Norm1()", xaxis, errorVSModelSize);
	}

}
