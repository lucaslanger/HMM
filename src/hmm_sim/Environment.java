package hmm_sim;

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
		FlowControl.outputData( "True_" + this.getDescription(), t);
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

	public int getDesiredHankelSize() {
		return this.desiredHankelSize;
	}


	
	
}
