package hmm_sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public class FlowControl {
	
	public static void main(String[] args){
		
		String workingFolder = "19-12-Labyrinth/";
		FlowControl.createFolder(workingFolder);
		int[] trajectorySizes = new int[]{25,50,100,200,500,1000,2000,4000,8000,16000};
		int repetitions = 200;
		
		//FlowControl.generateLabyrinthData(workingFolder, trajectorySizes, repetitions, 19, 12, 0, 200, .6, .4);
		
		int basisSize = 40;
		//FlowControl.readDataIntoModels(workingFolder, basisSize);
		
		testEngine a = new testEngine(workingFolder,"Models_Emperical_19_12_Toy_Labyrinth/", "Models_True_19_12_Toy_Labyrinth", 16000, basisSize, 2, 3);
		
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
		String inFolder = workingFolder + "Emperical_19_12_Toy_Labyrinth/";
		String outFolder = workingFolder + "Models_Emperical_19_12_Toy_Labyrinth/"; 
		FlowControl.createFolder(outFolder);
		
		String inFile = workingFolder + "True_19_12_Toy_Labyrinth";
		String TrueOut = workingFolder + "Models_True_19_12_Toy_Labyrinth";
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
	
	public static void generateLabyrinthData(String workingFolder, int[] trajectorySizes, int repetitions, int loop1, int loop2, double selfTransitionP , int hSize, double exit1P, double exit2P){
		rawHMM r = rawHMM.makeLabyrinth(workingFolder, loop1, loop2, selfTransitionP, hSize, exit1P, exit2P);
		r.generateData(trajectorySizes, repetitions);
	}
	
	public static void outputData(String s, double[][] data){
		ObjectOutputStream out;
		try {
			out = new ObjectOutputStream( new FileOutputStream( s ) );
			out.writeObject( data );
			out.close();
		} catch (IOException e) {
			System.out.println("Problem writing data");
			e.printStackTrace();
		}
		
		System.out.println("Done generating data for " + s);
		
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
