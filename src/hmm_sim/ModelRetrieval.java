package hmm_sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.Arrays;
import java.util.HashMap;

public class ModelRetrieval {
	
	
	private String fileNameOfTrueModel;
	private String workingFolder;
	private String fileNameOfEmpericalModels;
	private String pltFolder;
	private String[] fileNames;
	private int[] keySetSorted;
	private int basisSize;
	private int base;

	public ModelRetrieval(String workingFolder, String empModels,String fileNameOfTrueModel, int basisSize, int base){
		this.workingFolder = workingFolder;
		this.fileNameOfTrueModel = fileNameOfTrueModel;
		
		this.basisSize = basisSize;
		this.base = base;
		
		this.fileNameOfEmpericalModels = workingFolder + empModels;
		this.pltFolder = workingFolder + "Plotting_" + empModels + "/";
		testEngine.createFolder(this.pltFolder);
		
		File[] f = ModelRetrieval.getFiles(fileNameOfEmpericalModels);
		this.fileNames = new String[f.length];
		for (int i = 0; i < this.fileNames.length; i++) {
			this.fileNames[i] = f[i].getName();
		}
		this.makeKeySetSorted();
	}
	
	
	
	public HashMap<Integer, QueryEngine[]> getSpecificModelSizeQueryEngines(int numberOfTrajectoriesFromEachSize, int modelSize){
		HashMap<Integer, QueryEngine[]> dataSizeToModels = new HashMap<Integer, QueryEngine[]>();
		QueryEngine[] enginesTrajectorySize = new QueryEngine[numberOfTrajectoriesFromEachSize];
		QueryEngine q;

		try{
			for (String f: this.fileNames) {
				enginesTrajectorySize = new QueryEngine[numberOfTrajectoriesFromEachSize];	//Weird bug
				String file = this.fileNameOfEmpericalModels + f;
				int trajectoryLength = ModelRetrieval.getTrajectoryLengthFromFileName(file);
				ObjectInputStream ois = new ObjectInputStream( new FileInputStream(file) );
				HankelSVDModel h;
				
				for (int i = 0; i < numberOfTrajectoriesFromEachSize; i++) {
					h = (HankelSVDModel) ois.readObject();
					q = h.buildHankelBasedModel(this.basisSize, this.base, modelSize);
					enginesTrajectorySize[i] = q;
				}
				dataSizeToModels.put(trajectoryLength, enginesTrajectorySize);
				ois.close();
			}
			return dataSizeToModels;
		}
		catch(Exception e){
			System.out.println("Trouble making emperical models - fixed modelsize");
			e.printStackTrace();
			return null;
		}
	}
	
	private static File[] getFiles(String folder) {
		File dir = new File(folder);
		return dir.listFiles();
	}

	static void createFolder(String folder) {
		File dir = new File(folder);
		dir.mkdir();
		
	}
	
	private void makeKeySetSorted(){
		int[] kss = new int[this.fileNames.length];
		int i=0;
		for (String f: this.fileNames) {
			int trajectoryLength = ModelRetrieval.getTrajectoryLengthFromFileName(f);
			kss[i] = trajectoryLength;
			i++;
		}
		this.keySetSorted = kss;
	
		Arrays.sort(this.keySetSorted);
	}
	
	private static int getTrajectoryLengthFromFileName(String filename){
		return Integer.parseInt(filename.split(":")[1]);
	}

}
