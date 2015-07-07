package hmm_sim;

import java.util.Arrays;
import java.util.Random;

public abstract class Environment {
	
	
	private String desciption;
	private String empericalFolder;
	private String workingFolder;
	private String trueFile;
	private int desiredHankelSize;
	private int probabilityArraySize;
	private double[] trueProbabilities;

	public Environment(String workingFolder, String description, int desiredHankelSize){
		this.workingFolder = workingFolder;
		this.desciption = description;
		this.desiredHankelSize = desiredHankelSize;
		this.probabilityArraySize = 2*desiredHankelSize;
		this.empericalFolder = workingFolder + "Emperical_" + this.getDescription() + "/";
		this.trueFile = workingFolder + "True_" + this.getDescription();
	}
		
	public abstract double[] computeTrueProbabilities();
	
	public void initializeProbabilities(){
		this.trueProbabilities  = this.computeTrueProbabilities();
	}

	
	public double[] generateEmpericalProbabilities(int samples){
		double[] p = new double[this.probabilityArraySize];
		for (int i = 0; i < samples; i++) {
			int d = this.generateDuration();
			if (d < p.length){
				p[d] += 1;
			}
		}
		
		for (int i = 0; i < p.length; i++) {
			p[i] /= samples;
		}
		return p;
	}
	
	private int generateDuration(){ 						
		int l = this.trueProbabilities.length; 
		double[] cumulativeSum = new double[l];
		
		cumulativeSum[0] = this.trueProbabilities[0];
		for (int i = 1; i<l;i++){
			cumulativeSum[i] = cumulativeSum[i-1] + this.trueProbabilities[i];
		}

		Random random = new Random();
		double r = random.nextDouble();
		
		int index = Arrays.binarySearch(cumulativeSum, r);

		if (index >= 0){
			return index;
		}
		else{
			return -1*(index + 1);
		}
		
	}
	
	private double[][] makeHankel(double[] s){
		double[][] hankel = new double[this.desiredHankelSize][this.desiredHankelSize];
		for (int i = 0; i < this.desiredHankelSize; i++) {
			for (int j = 0; j < this.desiredHankelSize; j++) {
				hankel[i][j] = s[i+j];
			}
		}
	
		return hankel;
	}
	
	public void generateData(int[] trajectorySizes, int repetitions){
		FlowControl.createFolder(this.empericalFolder);
		this.printTrueProbabilities();
		for (int i = 0; i < trajectorySizes.length; i++) {
			this.printEmpericalTrials(trajectorySizes[i], repetitions);
		}
	}
	
	public void printTrueProbabilities(){
		//double[][] t = this.makeHankel(this.trueProbabilities);
		//FlowControl.outputData( this.getTrueFile(), t);
		FlowControl.outputData( this.getTrueFile(), new double[][]{ this.trueProbabilities });
	}
	
	public void printEmpericalTrials(int trajectoryLength, int repetitions){
		double[][] data = new double[repetitions][this.getDesiredHankelSize()*2];
		for (int i = 0; i < repetitions; i++) {
			data[i] = generateEmpericalProbabilities(trajectoryLength);
		}
		String f = this.empericalFolder + "Trajectory:" + Integer.toString(trajectoryLength);
		FlowControl.outputData(f, data );
	}


	public String getDescription() {
		return this.desciption;
	}
	
	public String getWorkingFolder(){
		return this.workingFolder;
	}
	
	public String getTrueFile(){
		return this.trueFile;
	}
	
	protected int getDesiredHankelSize() {
		return this.desiredHankelSize;
	}
	
	protected int getProbabilityArraySize() {
		return this.probabilityArraySize;
	}
	
	public double[] getTrueProbabilities(){
		return this.trueProbabilities;
	}
	
	
}
