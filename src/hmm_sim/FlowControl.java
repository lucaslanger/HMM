package hmm_sim;

import java.io.File;
import java.io.FileOutputStream;
import java.io.ObjectOutputStream;

public class FlowControl {
	
	public FlowControl(){
		
	}
	
	public static void readDataIntoModels(){
		int basisSize = 40;
		String inFolder = "Emperical_19_12_Toy_Labyrinth/";
		String outFolder = "Models_" + inFolder;
		File dir = new File(outFolder);
		dir.mkdir();
		
		File folder = new File(inFolder);
		File[] files = folder.listFiles();
		for (File file : files) {
			FlowControl.createModelsFromFile(inFolder, outFolder, file.getName(), "Models_" + file.getName(), basisSize);
		}
	}
	
	public static void createModelsFromFile(String inFolder, String outFolder, String fileIn, String fileOut, int basisSize){
		double[][] data = Environment.readData(inFolder + fileIn);
		HankelSVDModel[] modelsForFixedTrajectorySize = new HankelSVDModel[data.length];

		for (int i = 0; i < modelsForFixedTrajectorySize.length; i++) {
			modelsForFixedTrajectorySize[i] = new HankelSVDModel(data[i], basisSize);
		}
		
		FlowControl.outputModelsToFile(modelsForFixedTrajectorySize, outFolder, fileOut);
		
	}
	
	public static void outputModelsToFile(HankelSVDModel[] h, String outfolder, String fileOut){
		try{
			ObjectOutputStream oos = new ObjectOutputStream( new FileOutputStream(outfolder + fileOut) );
			System.out.println(outfolder + fileOut);
			for (int i = 0; i < h.length; i++) {
				oos.writeObject(h[i].getProbabilities());
				oos.writeObject(h[i].getBasisSize());
				oos.writeObject(h[i].getSvd());
			}
			oos.close();
		}
		catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void generateLabyrinthData(){
		int[] trajectorySizes = new int[]{25,50,100,200,500,1000,2000,4000,8000,16000};
		int repetitions = 100;
		rawHMM r = rawHMM.makeLabyrinth(19, 12, 0, 200);
		r.generateData(trajectorySizes, repetitions);
	}
	
	public static void createFolder(String folder){
		File dir = new File(folder);
		dir.mkdir();
	}

}
