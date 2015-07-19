package hmm_sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

import Jama.Matrix;

public class FlowControl {
	
	public static void main(String[] args){
		int[] trajectorySizes = new int[]{25,50,100,200,500,1000,2000,4000,8000,16000,32000,64000,128000,256000};
		int dataSizeForFixedPlots = 256000;
		int base = 2; // Haven't tested for bases other than 2 ... no guarantees
	
		//FlowControl.testLabyrinths(trajectorySizes, dataSizeForFixedPlots, base);
		//FlowControl.testLoops(trajectorySizes, dataSizeForFixedPlots, base);
		FlowControl.computeKeySearchStuff(trajectorySizes, dataSizeForFixedPlots, base);
	}
	
	public static void testLabyrinths(int[] trajectorySizes, int dataSizeForFixedPlots, int base){
		
		int repetitions = 5;
		int stretchFactor = 10;
		int hSize = 500;
		int basisSize = 300;
		int[] modelSizes = new int[]{};

		String workingFolder = "testLargeLabyrinth/";
	
		System.out.println("Generating data:");
		System.out.println("");
		FlowControl.createFolder(workingFolder);
		LabyrinthGraph l = LabyrinthGraph.pacMan(workingFolder, hSize, stretchFactor, 5);
		//LabyrinthGraph l = LabyrinthGraph.testLabyrinth(workingFolder, hSize, stretchFactor);
		l.generateData(trajectorySizes, repetitions, false);
		
		System.out.println("");
		
		System.out.println("Reading data into models");
		FlowControl.readDataIntoModels(workingFolder, basisSize);
		System.out.println("Done loading models");
		System.out.println("");
		
		testEngine a = new testEngine(workingFolder,"Models_Emperical_" + workingFolder, "Models_True_" + workingFolder, dataSizeForFixedPlots, basisSize, base, modelSizes, 50 ,1, true );
	}
	
	public static void testLoops(int[] trajectorySizes, int dataSizeForFixedPlots, int base){
		int repetitions = 50;

		int loop1 = 64;
		int loop2 = 16;
		int hSize = 450;
		int basisSize = 150;
		int[] modelSizes = new int[]{};
		//Bug of having all errors exactly the same seems to occur when taken model size is really large e.g 50 was tried
		
		String workingFolder = Integer.toString(loop1) + "_" + Integer.toString(loop2) + "_Toy_Labyrinth/";
		
		System.out.println("Generating data:");
		System.out.println("");
		FlowControl.createFolder(workingFolder);
		rawHMM r = rawHMM.makeLabyrinth(workingFolder, loop1, loop2, 0.10, hSize, .5, .5);
		r.generateData(trajectorySizes, repetitions, false);
		
		System.out.println("");
		
		System.out.println("Reading data into models");
		FlowControl.readDataIntoModels(workingFolder, basisSize);
		System.out.println("Done loading models");
		System.out.println("");
		
		testEngine a = new testEngine(workingFolder,"Models_Emperical_" + workingFolder, "Models_True_" + workingFolder, dataSizeForFixedPlots , basisSize, base, modelSizes, 30, 2 , true);
	}
	
	public static void computeKeySearchStuff(int[] trajectorySizes, int dataSizeForFixedPlots, int base){
		int repetitions = 5;
		int stretchFactor = 10;
		int hSize = 500;
		int basisSize = 300;
		int key = 5;
				
		int samples = 100;
		String workingFolder = "keySearchPacMan/";
		String empModels = "Models_Emperical_" + workingFolder;
		String pltFolder = workingFolder + "Plotting_" + empModels + "/";
		
		FlowControl.createFolder(workingFolder);

		LabyrinthGraph l = LabyrinthGraph.pacMan(workingFolder, hSize, stretchFactor, key);
		l.generateData(trajectorySizes, repetitions, false);
		
		FlowControl.readDataIntoModels(workingFolder, basisSize);

		double[][] modelSizes = new double[][]{{3, 5, 10, 20, 30, 50, 70, 90}};
		double[][] maxKsToTest = new double[][]{{5, 10, 20, 40, 80, 160, 320}};
		double[][] errorVSModelSize = new double[maxKsToTest[0].length][modelSizes[0].length];
		double[][] xaxis = new double[maxKsToTest[0].length][modelSizes[0].length];

		int[] shortestPaths = l.shortestPathsFromKey();

		System.out.println("AVERAGING");
		System.out.println("HACK ON PSEUDOINVERSE");
			
		for (int j = 0; j < maxKsToTest[0].length; j++) {
			int k = (int) maxKsToTest[0][j]; 
			int[][] durationDistancePairs = l.createObservationDistanceSamples(shortestPaths, k, samples);
			double[][] trueDistanceAhead = l.dynamicallyDetermineTrueDistanceKAhead(shortestPaths, k);
			
			for (int i = 0; i < modelSizes[0].length; i++) {		
				int m = (int) modelSizes[0][i];
				
				//testEngine a = new testEngine(workingFolder, empModels, "Models_True_" + workingFolder, dataSizeForFixedPlots , basisSize, base, new int[]{}, m, 1 , false);
				//QueryEngine learnedModel = a.fixedModelQE.get(dataSizeForFixedPlots)[0];

				ModelRetrieval mr = new ModelRetrieval(workingFolder, empModels, "Models_True_" + workingFolder, basisSize, base);
				QueryEngine learnedModel = mr.getSpecificModelSizeQueryEngines(1, m).get(dataSizeForFixedPlots)[0];
				
				Matrix[] alphaKStates = learnedModel.getAllKStateQueries(k, base);
				
				Matrix Atheta = l.getAlphaFromSampledData(durationDistancePairs, alphaKStates);
	
				double e = l.performanceDistanceErrorComputations(Atheta, trueDistanceAhead, durationDistancePairs);
				errorVSModelSize[j][i] = e;
				xaxis[j][i] = i;
			}
			
		}
		Matrix pe = new Matrix(errorVSModelSize);
		pe.print(5, 5);
		
		System.out.println("Writing data to: " + pltFolder + "KeyFindingError");
		OutputData.outputData(pltFolder + "KeyFindingError", "modelSize", "error", xaxis, errorVSModelSize);
	}
	
	public FlowControl(){
		
	}
	
	public static void createFolder(String folder){
		File dir = new File(folder);
		dir.mkdir();
	}
	
	public static File[] getFiles(String folder){
		File dir = new File(folder);
		return dir.listFiles();
	}
	
	public static void readDataIntoModels(String workingFolder, int basisSize){
		String inFolder = workingFolder + "Emperical_" + workingFolder;
		String outFolder = workingFolder + "Models_Emperical_" + workingFolder; 
		FlowControl.createFolder(outFolder);
		
		String inFile = workingFolder + "True_" + workingFolder;
		String TrueOut = workingFolder + "Models_True_" + workingFolder;
		FlowControl.createModelsFromFile("", "", inFile, TrueOut, basisSize);
		
		File[] files = FlowControl.getFiles(inFolder);
		for (File file : files) {
			FlowControl.createModelsFromFile(inFolder, outFolder, file.getName(), "Models_" + file.getName(), basisSize);
		}
	}	
	
	
	public static void createModelsFromFile(String inFolder, String outFolder, String fileIn, String fileOut, int basisSize){
		double[][] data = FlowControl.readData(inFolder + fileIn);
		HankelSVDModel[] modelsForFixedTrajectorySize = new HankelSVDModel[data.length];

		for (int i = 0; i < modelsForFixedTrajectorySize.length; i++) {
			modelsForFixedTrajectorySize[i] = new HankelSVDModel(data[i], basisSize);
		}
		
		FlowControl.outputModelsToFile(modelsForFixedTrajectorySize, outFolder, fileOut);
		
	}
	
	private static void outputModelsToFile(HankelSVDModel[] h, String outfolder, String fileOut){
		try{
			ObjectOutputStream oos = new ObjectOutputStream( new FileOutputStream(outfolder + fileOut) );
			for (int i = 0; i < h.length; i++) {
				oos.writeObject(h[i]);
			}
			oos.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
		
	public static void outputData(String s, double[][] data, boolean verbose){
		ObjectOutputStream out;
		try {
			out = new ObjectOutputStream( new FileOutputStream( s ) );
			out.writeObject( data );
			out.close();
		} catch (IOException e) {
			System.out.println("Problem writing data");
			e.printStackTrace();
		}
		
		if(verbose){
			System.out.println("Done generating data for " + s);
		}
	}
	
	public static double[][] readData(String filename){
		double[][] data;
		try{
			ObjectInputStream ois = new ObjectInputStream( new FileInputStream(filename) );
			data = (double[][]) ois.readObject();
			ois.close();
			return data;
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
	}


}
