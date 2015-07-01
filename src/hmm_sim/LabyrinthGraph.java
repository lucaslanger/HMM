package hmm_sim;

import Jama.Matrix;


public class LabyrinthGraph extends Environment{
	
	private int[][] graph;
	private int[][] edges;
	private double[][] transitions; 
	private double[] prior;
	private double[] transitionsToDoor;
		
	public LabyrinthGraph(String workingFolder, String description, int desiredHankelSize, int[][] graph, int[][] edges, double[][] transitions, double[] transitionsToDoor, double[] prior){
		super(workingFolder, description, desiredHankelSize);
		this.graph = graph;
		this.edges = edges;
		this.transitions = transitions;
		this.transitionsToDoor = transitionsToDoor;
		this.prior = prior;
		this.computeTrueProbabilities();
		
	}
	
	public double[] computeTrueProbabilities(){
		int iterations = 100;
		
		int c=0;
		
		double[][] nodeToDuration = new double[this.graph.length][super.getProbabilityArraySize()];
		
		while (c<iterations) {
			for (int i = 0; i < this.graph.length; i++) {
				computeLocalDurations(i, nodeToDuration);
			}
		}
		
		Matrix n = new Matrix(nodeToDuration);
		n.print(5, 5);
		n.transpose().print(5, 5);
		
		Matrix p = new Matrix(prior, 0) ;
		Matrix r = n.times( p );
		return r.getArrayCopy()[0];
	}
	
	
	public void computeLocalDurations(int node, double[][] nodeToDuration){
		
		double[] rDurationDistribution = new double[super.getProbabilityArraySize()];

		for (int i = 0; i < this.graph[node].length; i++) {
			int neighbor = this.graph[node][i];
			if (neighbor == -1) {
				rDurationDistribution[0] += this.transitionsToDoor[node];
				continue;
			} 
			
			else{ 
				int transitionDuration = this.edges[node][i];
				double transitionProbability = this.transitions[node][i];
				double[] contribution = LabyrinthGraph.shiftContribution(transitionDuration, transitionProbability, nodeToDuration[neighbor] );
				rDurationDistribution = LabyrinthGraph.addDoubleArrays(rDurationDistribution, contribution);
			}
		}
		
		nodeToDuration[node] = rDurationDistribution;
		
	}

	private static double[] shiftContribution(int transitionDuration, double transitionProbability, double[] durationProbabilities) {
		double[] rDurations = new double[durationProbabilities.length];
		for (int i = 0; i < durationProbabilities.length - transitionDuration; i++) {
			rDurations[transitionDuration + i] = durationProbabilities[i] * transitionProbability;
		}
		return rDurations;				
	}
	
	private static double[] addDoubleArrays(double[] d1, double[] d2){
		double[] d3 = new double[d2.length];
		for (int i = 0; i < d2.length; i++) {
			d3[i] = d1[i] + d2[i];
		}
		return d3;
	}
}


/*
public double[] computeTrueProbabilities(){
	double[] probabilities = new double[super.getProbabilityArraySize()];
	
	HashMap<Integer, double[]> preComputed = new HashMap<Integer, double[]>();
	
	for (int i = 0; i < this.prior.length; i++) {
		double[] t = new double[probabilities.length];
		
		depthLimitedSearch(i, this.prior[i], t, 0, super.getProbabilityArraySize(), preComputed);
		for (int j = 0; j < t.length; j++) {
			probabilities[j] += t[j];
		}
		
	}
	return probabilities;
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
*/
