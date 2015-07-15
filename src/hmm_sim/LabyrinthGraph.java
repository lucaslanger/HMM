package hmm_sim;

import java.awt.AlphaComposite;
import java.util.Arrays;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Random;

import Jama.Matrix;


public class LabyrinthGraph extends Environment{
	
	private int[][] graph;
	private int[][] edges;
	private double[][] transitions; 
	private double[] prior;
		
	public LabyrinthGraph(String workingFolder, int desiredHankelSize, int[][] graph, int[][] edges, double[][] transitions, double[] prior){
		super(workingFolder, workingFolder, desiredHankelSize);
		this.graph = graph;
		this.edges = edges;
		this.transitions = transitions;
		this.prior = prior;
		super.initializeProbabilities();
		
	}
	
	public double[] computeTrueProbabilities(){
		int iterations = 100;
		
		int c=0;
		
		double[][] nodeToDuration = new double[this.graph.length][super.getProbabilityArraySize()];
		
		while (c<iterations) {
			for (int i = 1; i < this.graph.length; i++) {
				computeLocalDurations(i, nodeToDuration);
			}
			c++;
		}
		
		Matrix n = new Matrix(nodeToDuration).transpose();
		
		Matrix p = new Matrix( new double[][]{ prior } ).transpose() ;
		Matrix r = n.times( p );
		
		System.out.println("Probabilities:");
		r.transpose().print(5, 5);
		System.out.println("Sum");
		System.out.println(LabyrinthGraph.sumArray(r.transpose().getArrayCopy()[0]));
		System.out.println();
		return r.transpose().getArrayCopy()[0];
	}
	
	
	public void computeLocalDurations(int node, double[][] nodeToDuration){
		
		double[] rDurationDistribution = new double[super.getProbabilityArraySize()];

		for (int i = 0; i < this.graph[node].length; i++) {
			int neighbor = this.graph[node][i];
			if (neighbor == 0) {
				rDurationDistribution[0] += this.transitions[node][0];
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
	
	private static double sumArray(double[] da){
		double s = 0;
		for(double d: da){
			s += d;
		}
		return s;
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
	
	public static LabyrinthGraph pacMan(String workingFolder, int desiredHankelSize, int stretchFactor ){
			
		 int[][] graph = new int[][]{   
				 {},
				 {0,2},
				 {0,3},
				 {4,13},
				 {5,13},
				 {6,15},
				 {0,7},
				 {0,8,16},
				 {9,16},
				 {0,10},
				 {0,11},
				 {12,14},
				 {0,13},
				 {0,3,4},
				 {11,15,16},
				 {5,14},
				 {7,8}
		 };
		 int[][] edges = new int[][]{
				 {},
				 {0,1},
				 {0,1},
				 {2,3},
				 {1,2},
				 {1,2},
				 {0,1},
				 {0,2,1},
				 {1,1},
				 {0,1},
				 {0,3},
				 {2,2},
				 {0,1},
				 {0,1,2},
				 {1,1,2},
				 {2,3},
				 {2,1}
		 };
		 
		 LabyrinthGraph.stretchEdges(edges, stretchFactor);

		 double[][] transitions = new double[edges.length][];
		 for (int i = 0; i < edges.length; i++) {
			transitions[i] = LabyrinthGraph.equalVector(edges[i].length);
								//LabyrinthGraph.randomVector(edges[i].length);
		}
		 
		 
		 double[] prior = new double[edges.length];
		 prior[5] = 1;

		 LabyrinthGraph l = new LabyrinthGraph(workingFolder, desiredHankelSize, graph, edges, transitions, prior);
		 return l;
	}
	
	private static double[] equalVector(int length) {
		double[] d = new double[length];
		for (int i = 0; i < d.length; i++) {
			d[i] = 1.0/ d.length;
		}
		return d;
	}

	private static double[] randomVector(int length) {
		double[] d = new double[length];
		Random r = new Random();
		double sum = 0;
		for (int i = 0; i < d.length; i++) {
			d[i] = r.nextDouble();
			sum += d[i];
		}
		for (int i = 0; i < d.length; i++) {
			d[i] /= sum;
		}
		return d;
	}

	public static LabyrinthGraph testLabyrinth(String workingFolder, int desiredHankelSize, int stretchFactor){
		
		 int[][] graph = new int[][]{   
				 {},
				 {0, 2},
				 {1, 3},
				 {2, 4},
				 {0, 3}
		 };
		 int[][] edges = new int[][]{
				 {},
				 {0, 1},
				 {2, 2},
				 {1, 3},
				 {0, 2}
		 };
		 
		 LabyrinthGraph.stretchEdges(edges, stretchFactor);

		 double[][] transitions = new double[][]{
				 {},
				 {0.5, 0.5},
				 {0.5, 0.5},
				 {0.5, 0.5},
				 {0.5, 0.5}
		 };
		 
		 
		 double[] prior = new double[]{0,1,0,0,0};

		 LabyrinthGraph l = new LabyrinthGraph(workingFolder, desiredHankelSize, graph, edges, transitions, prior);
		 return l;
		 
	}
	
	private static void stretchEdges(int[][] edges, int stretchFactor) {
		for (int i = 0; i < edges.length; i++) {
			for (int j = 0; j < edges[i].length; j++) {
				edges[i][j] = edges[i][j]*stretchFactor;
			}
		}
	}

	public static void main(String[] args){
		String wf = "test";
		int d = 100;
		LabyrinthGraph.testLabyrinth(wf,  d, 1);
	}
	
	
	public HashMap<Integer, Integer> shortestPathsFromKey(){
		int start = 5;
		
		HashMap<Integer, Integer> paths = new HashMap<Integer, Integer>();
		HashMap<Integer, DijkstraNode> accessIntoHeap = new HashMap<Integer, DijkstraNode>();
		PriorityQueue<DijkstraNode> pq = new PriorityQueue<DijkstraNode>();
		
		DijkstraNode startNode = new DijkstraNode(start, 0);
		pq.add(startNode);
		
		while(pq.isEmpty() == false){
			DijkstraNode n = pq.remove();
			int lengthToN = n.getLengthToNode();
			paths.put(n.getId(), lengthToN) ;
			for (int i = 0; i < graph[n.getId()].length; i++) {
				int outEdge = graph[n.getId()][i];
				DijkstraNode d;
				if ( accessIntoHeap.containsKey(outEdge) ){
					int l = lengthToN + edges[n.getId()][i];
					d = accessIntoHeap.get(outEdge);
					if ( d.getLengthToNode() < l ){
						d.setLengthToNode(l);
					}
				}
				else{
					d = new DijkstraNode(outEdge, lengthToN + edges[n.getId()][i]);
				}
				pq.add(d);
			}
		}

		return paths;
	}
	
	public Matrix getAlphaFromSampledData(int[][] samples, Matrix[] alphaKStates){
		Matrix distances = new Matrix(samples)[0];
		Matrix durations = new Matrix(samples)[1];
		
		double[][] a = new double[samples.length][alphaKStates[0].getArray().length];
		for (int i = 0; i < samples.length; i++) {
			int k = (int) durations.get(i, 0);
			double[] aK = alphaKStates[k].getArrayCopy()[0];
			a[i] = aK;
		}
		Matrix A = new Matrix(a);
		Matrix AT = A.transpose();
		Matrix inverseForRegression = AT.times(A).inverse();
		
		Matrix theta = inverseForRegression.times(distances);
		return theta;
	}
	
	public performanceDistanceErrorComputations(Matrix theta, Matrix[] alphaKStates, double[][] trueDistanceKAhead, int[][] samples){
		double[] trueAverageDistance = new double[samples.length];
		for (int i = 0; i < trueAverageDistance.length; i++) {
			trueAverageDistance[i] = trueDistanceKAhead[samples[i][0]];
		}
		A.times(theta);
	}
	
	public int[][] createObservationDistanceSamples(HashMap<Integer, Integer> shortestPaths,int maxObservation, int samples, double[] initialStates){
		int[][] s = new int[samples][2];
		Random random = new Random();
		for (int i = 0; i < samples; i++) {
			int r = random.nextInt(maxObservation); 
			int d = sampleDistance(shortestPaths, r, random ,initialStates);
			s[i] = new int[]{r,d};
		}
		return s;
	}

	private int sampleDistance(HashMap<Integer, Integer> shortestPaths, int r, Random random, double[] initialStates) {
		while(true){
			int state = sampleState(initialStates, random);
			int d = r;
			while(true){
				int desiredNextState = graph[state][sampleState(transitions[state], random)];
				if (desiredNextState == 0){
					break;
				}
				else if( edges[state][desiredNextState] > d){
					return shortestPaths.get(state) + d;
				}
				else{
					d = d - edges[state][desiredNextState]; 
					state = desiredNextState;
				}
			}
		}
	}
	
	public int sampleState(double[] d, Random random){
		double[] c = cumulativeSum(d);
		double r = random.nextDouble();
		
		int index = Arrays.binarySearch(c, r);
		return index;
	}
	
	public double[] cumulativeSum( double[] d){
		double[] r = new double[d.length];
		r[0] = d[0];
		for (int i = 1; i < r.length; i++) {
			r[i] = r[i-1] + d[i]; 
		}
		return r;
	}
	
	private double[][] dynamicallyDetermineTrueDistanceKAhead(HashMap<Integer, Integer> shortestPaths, int maxK){
		double[][] distanceStorage = new double[maxK][this.graph.length];		
		
		for (int i = 0; i < maxK; i++) {
			for (int j = 0; j < distanceStorage.length; j++) {
				if (i==0){
					distanceStorage[i][j] = shortestPaths.get(j);
				}
				else{
					for (int j2 = 0; j2 < graph[j].length; j2++) {
						int n = graph[j][j2];
						int d = i - edges[j][j2];
						if(d >= 0){
							distanceStorage[i][j] += distanceStorage[d][n]*transitions[j][j2];
						}
						else{
							distanceStorage[i][j] += (distanceStorage[0][n] + -1*d)*transitions[j][j2];
						}
					}
				}
				
			}
		}
		return distanceStorage;
	}

}
