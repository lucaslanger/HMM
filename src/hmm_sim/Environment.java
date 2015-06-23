package hmm_sim;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public abstract class Environment {
	
	
	private String desciption;
	public String empericalFolder; 
	private int desiredHankelSize;

	public Environment(String description, int desiredHankelSize){
		this.desciption = description;
		this.desiredHankelSize = desiredHankelSize;
		this.empericalFolder = "Emperical_" + this.getDescription() + "/";
	}
		
	public abstract double[] generateEmpericalProbabilities(int samples);
	
	public abstract double[][] generateTrueProbabilities();
	
	public void generateData(int[] trajectorySizes, int repetitions){
		FlowControl.createFolder(this.empericalFolder);
		this.printTrueProbabilities();
		for (int i = 0; i < trajectorySizes.length; i++) {
			this.printEmpericalTrials(trajectorySizes[i], repetitions);
		}
	}
	
	public void printTrueProbabilities(){
		double[][] t = this.generateTrueProbabilities();
		this.outputData( "TrueModel_", t);
	}
	
	public void printEmpericalTrials(int trajectoryLength, int repetitions){
		double[][] data = new double[repetitions][this.getDesiredHankelSize()*2];
		for (int i = 0; i < repetitions; i++) {
			data[i] = generateEmpericalProbabilities(trajectoryLength);
		}
		this.outputData(this.empericalFolder + "Trajectory:" + Integer.toString(trajectoryLength), data );
	}
	
	public void outputData(String s, double[][] data){
		ObjectOutputStream out;
		try {
			out = new ObjectOutputStream( new FileOutputStream( s ) );
			out.writeObject( data );
			out.close();
		} catch (IOException e) {
			System.out.println("Problem writing data");
			e.printStackTrace();
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


	public String getDescription() {
		return this.desciption;
	}

	public int getDesiredHankelSize() {
		return this.desiredHankelSize;
	}


	
	
}
