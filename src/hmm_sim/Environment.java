package hmm_sim;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;

public abstract class Environment {
	
	private String description;
	private int desiredHankelSize;
	
	public abstract double[] generateEmpericalProbabilities(int samples);
	
	public abstract double[][] generateTrueProbabilities();
	
	public void printTrueProbabilities(){
		double[][] t = this.generateTrueProbabilities();
		this.outputData( generateTrueProbabilities() );
	}
	
	public void printEmpericalTrials(int numTrials, int samples){
		double[][] data = new double[numTrials][samples];
		for (int i = 0; i < numTrials; i++) {
			data[i] = generateEmpericalProbabilities(samples);
		}
		this.outputData(data);
	}
	
	public void outputData(double[][] data){
		ObjectOutputStream out;
		try {
			out = new ObjectOutputStream( new FileOutputStream(this.description) );
			out.writeObject(data);
			out.flush();
			out.close();
		} catch (IOException e) {
			System.out.println("Problem writing data");
			e.printStackTrace();
		}
		
	}
	
	public static double[][] readData(String filename) throws FileNotFoundException, IOException, ClassNotFoundException{
		double[][] data;
		
		ObjectInputStream ois = new ObjectInputStream( new FileInputStream(filename) );
		data = (double[][]) ois.readObject();
	
		return data;
		
	}
	
	public abstract String getDescription();
	
	public abstract int getDesiredHankelSize();
	
	
}
