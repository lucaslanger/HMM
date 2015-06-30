package hmm_sim;

import java.util.HashMap;


public class LabyrinthGraph extends Environment{
	
	private int[][] graph;
	private int[][] edges;
	private double[][] transitions; 
	private double[] prior;
	
	private double[] trueProbabilities;
	
	public LabyrinthGraph(String workingFolder, String description, int desiredHankelSize, int[][] graph, int[][] edges, double[][] transitions, double[] prior){
		super(workingFolder, description, desiredHankelSize);
		this.graph = graph;
		this.edges = edges;
		this.transitions = transitions;
		this.prior = prior;
		this.computeTrueProbabilities();
	}
	
	public void computeTrueProbabilities(){
		int pSize  = super.getDesiredHankelSize()*2;
		double[] probabilities = new double[pSize];
		
		HashMap<Integer, double[]> preComputed = new HashMap<Integer, double[]>();
		
		for (int i = 0; i < this.prior.length; i++) {
			double[] t = new double[probabilities.length];
			
			depthLimitedSearch(i, this.prior[i], t, 0, pSize, preComputed);
			for (int j = 0; j < t.length; j++) {
				probabilities[j] += t[j];
			}
			
		}
		this.trueProbabilities = probabilities;
	}
	
	private void depthLimitedSearch(int currentNode, double p, double[] terminations, int currentDepth, int depthLimit, HashMap<Integer, double[]> alreadyComputed){
		if (p==0 || currentDepth >= depthLimit){
			//DoNothing	
		}
		else if(currentNode == -1){
			terminations[currentDepth] += p;
		}
		for (int i = 0; i < this.graph[currentNode].length; i++) {
			int n = this.graph[currentNode][i];
			depthLimitedSearch(n, p*this.transitions[currentNode][n], terminations, currentDepth + this.edges[currentNode][n], depthLimit, alreadyComputed);
		}
	}


	private double[][] makeHankel( double[] probabilities){
		double[][] hankel = new double[probabilities.length/2][probabilities.length/2];
		for (int i = 0; i < probabilities.length; i++) {
			for (int j = 0; j < hankel.length; j++) {
				hankel[i][j] = probabilities[i+j];
			}
		}
		return hankel;
	}

	@Override
	public double[] generateTrueProbabilities() {
		return this.trueProbabilities;
	}
	
}
