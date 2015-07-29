package hmm_sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.util.ArrayList;
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
	
	public ModelEnginePair getAllSizesQueryEngines(int numberOfTrajectoriesFromEachSize, int numberOfModels, int[] modelSizes){
		
		HashMap<Integer, HankelSVDModel[]> hsvdModels = new HashMap<Integer, HankelSVDModel[]>();
		HashMap<Integer, QueryEngine[][]> dataSizeToModels = new HashMap<Integer, QueryEngine[][]>();
		QueryEngine[][] enginesModelSizeTrajectorySize = new QueryEngine[numberOfModels][numberOfTrajectoriesFromEachSize];
		QueryEngine q;

		try{
			for (String f: this.fileNames) {
				enginesModelSizeTrajectorySize = new QueryEngine[numberOfModels][numberOfTrajectoriesFromEachSize];	//Weird bug
				String file = this.fileNameOfEmpericalModels + f;
				int trajectoryLength = ModelRetrieval.getTrajectoryLengthFromFileName(file);
				ObjectInputStream ois = new ObjectInputStream( new FileInputStream(file) );
				
				HankelSVDModel[] ha = new HankelSVDModel[numberOfTrajectoriesFromEachSize];
				for (int i = 0; i < numberOfTrajectoriesFromEachSize; i++) {
					HankelSVDModel h = new HankelSVDModel();
					h = (HankelSVDModel) ois.readObject();
					for (int j = 0; j < numberOfModels; j++) {
						q = h.buildHankelBasedModel(this.base, modelSizes[j]);
						enginesModelSizeTrajectorySize[j][i] = q;
					}
					ha[i] = h;
				}
				dataSizeToModels.put(trajectoryLength, enginesModelSizeTrajectorySize);
				hsvdModels.put(trajectoryLength, ha);
				ois.close();
			}
			
			return new ModelEnginePair(dataSizeToModels, hsvdModels);
			
		}
		catch(Exception e){
			System.out.println("Trouble making emperical models");
			e.printStackTrace();
			return null;
		}
		
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
					q = h.buildHankelBasedModel(this.base, modelSize);
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
	
	public HankelSVDModel readTrueModel(String f){
		f = this.workingFolder + f;
		HankelSVDModel t;
		//Search for file with True model
		ObjectInputStream ois;
		try {
			ois = new ObjectInputStream(new FileInputStream(f));
			t = (HankelSVDModel) ois.readObject();
			ois.close();
			return t;
		} catch (IOException e) {
			System.out.println("Problem fetching true Model");
			e.printStackTrace();
			return null;
		} catch (ClassNotFoundException e) {
			System.out.println("Couldn't find the file");
			e.printStackTrace();
			return null;
		}
		
	}
	
	public double[] checkEngine(QueryEngine q, HankelSVDModel trueModel, String id, int topCount) {
		System.out.println("Inspecting the following engine: " + id);
		System.out.println();
		
		double[] error = new double[q.getMaxExponent()+1];
		for (int i = 0; i <= q.getMaxExponent(); i++) {
			int pow = (int) Math.pow(this.base,i);
		
			double[] r = new double[trueModel.getProbabilities().length];
			for (int i1 = 0; i1 < r.length; i1++) {
				r[i1] = Math.abs(trueModel.getProbabilities()[i1] - q.probabilityQuery(i1, pow, this.base, true));
			}
			System.out.print ("Qerrors maxpower: " + Integer.toString(pow) + ":  ");
			System.out.print( Arrays.toString(getTopErrorIndices(r, topCount)[0]) );
			System.out.print(", ");
			System.out.println( Arrays.toString(getTopErrorIndices(r, topCount)[1]) );
			q.debugProbabilityQuery(160, pow, this.base, true);

			double s = ModelRetrieval.sumArray(r);
			System.out.println(s);
			error[i] = s; 
			
		}
		System.out.println();
		return error;
	}
	
	private static int getTrajectoryLengthFromFileName(String filename){
		return Integer.parseInt(filename.split(":")[1]);
	}
	
	public static double[][] getTopErrorIndices(double[] r, int topcount){
		double[] d = new double[r.length];
		for (int i = 0; i <r.length; i++) {
			d[i] = r[i];
		}
		HashMap<Double, ArrayList<Integer>> val_to_index = new HashMap<Double, ArrayList<Integer>>();
		for (int i = 0; i < d.length; i++) {
			ArrayList<Integer> v = val_to_index.get(d[i]);
			if (v == null){
				v = new ArrayList<Integer>();
			}
			v.add(i);
			val_to_index.put(d[i], v);
		}
		Arrays.sort(d);
		double[] reverse_d = new double[d.length];
		for (int i = 0; i < d.length; i++) {
			reverse_d[d.length-1-i] = d[i];
		}
		d = reverse_d;
	
		double[] top = new double[topcount];
		double[] topvalues = new double[topcount];
		
		int i = 0;
		int t = 0;
		while(t < top.length){
			while(t-i<val_to_index.get(d[i]).size() && t < top.length){
				top[t] = val_to_index.get(d[i]).get(t-i);
				topvalues[t] = d[i];
				t++;
			}
			i++;
		}
		return new double[][]{ top, topvalues };
		
		
	}
	
	private static double sumArray(double[] da){
		double s = 0;
		for(double d: da){
			s += d;
		}
		return s;
	}
	
	public String getFileNameOfTrueModel() {
		return fileNameOfTrueModel;
	}



	public String getWorkingFolder() {
		return workingFolder;
	}



	public String getFileNameOfEmpericalModels() {
		return fileNameOfEmpericalModels;
	}



	public String getPltFolder() {
		return pltFolder;
	}



	public String[] getFileNames() {
		return fileNames;
	}



	public int[] getKeySetSorted() {
		return keySetSorted;
	}



	public int getBasisSize() {
		return basisSize;
	}



	public int getBase() {
		return base;
	}

}
