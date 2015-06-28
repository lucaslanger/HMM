package hmm_sim;

public abstract class Environment {
	
	
	private String desciption;
	private String empericalFolder;
	private String workingFolder;
	private String trueFile;
	private int desiredHankelSize;

	public Environment(String workingFolder, String description, int desiredHankelSize){
		this.workingFolder = workingFolder;
		this.desciption = description;
		this.desiredHankelSize = desiredHankelSize;
		this.empericalFolder = workingFolder + "Emperical_" + this.getDescription() + "/";
		this.trueFile = workingFolder + "True_" + this.getDescription();
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
		FlowControl.outputData( this.getTrueFile(), t);
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

	public int getDesiredHankelSize() {
		return this.desiredHankelSize;
	}
	


	
	
}
