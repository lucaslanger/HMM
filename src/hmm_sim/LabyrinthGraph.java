package hmm_sim;

import java.awt.AlphaComposite;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Random;

import Jama.Matrix;


public class LabyrinthGraph extends Environment{
	
	private int[][] graph;
	private int[][] edges;
	private double[][] transitions; 
	private double[] prior;
	private int key;
	private HashMap<Integer, ArrayList<Integer>> incomingEdges;
		
	public LabyrinthGraph(String workingFolder, int desiredHankelSize, int[][] graph, int[][] edges, double[][] transitions, double[] prior){
		super(workingFolder, workingFolder, desiredHankelSize);
		this.graph = graph;
		this.edges = edges;
		this.transitions = transitions;
		this.prior = prior;
		super.initializeProbabilities();
		
	}
	
	public LabyrinthGraph(String workingFolder, int desiredHankelSize, int[][] graph, int[][] edges, double[][] transitions, double[] prior, int key){
		super(workingFolder, workingFolder, desiredHankelSize);
		this.graph = graph;
		this.edges = edges;
		this.transitions = transitions;
		this.prior = prior;
		this.key = key;
		this.buildIncomingEdges();
		super.initializeProbabilities();
		
	}
	
	public void buildIncomingEdges(){
		HashMap<Integer, ArrayList<Integer>> incomingEdges = new HashMap<Integer, ArrayList<Integer>>();
		for (int i = 0; i < graph.length; i++) {
			for (int j = 0; j < graph[i].length; j++) {
				if (incomingEdges.get( graph[i][j] ) != null){
					incomingEdges.get( graph[i][j] ).add(i);
				}
				else{
					ArrayList<Integer> a = new ArrayList<Integer>();
					a.add(i);
					incomingEdges.put( graph[i][j], a );
				}
			}
		}
		this.incomingEdges = incomingEdges;
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
	
	public static LabyrinthGraph pacMan(String workingFolder, int desiredHankelSize, int stretchFactor, int key){
			
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
		 prior[3] = 1;

		 LabyrinthGraph l = new LabyrinthGraph(workingFolder, desiredHankelSize, graph, edges, transitions, prior, key);
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
	
	
	public int[] shortestPathsFromKey(){
		int key = this.key;
		
		HashMap<Integer, Integer> paths = new HashMap<Integer, Integer>();
		HashMap<Integer, DijkstraNode> accessIntoHeap = new HashMap<Integer, DijkstraNode>();
		PriorityQueue<DijkstraNode> pq = new PriorityQueue<DijkstraNode>( 20, new DijkstraNodeComparator() );
		
		DijkstraNode startNode = new DijkstraNode(key, 0);
		pq.add(startNode);
		
		while(pq.isEmpty() == false){
			DijkstraNode n = pq.remove();
			if (paths.containsKey(n.getId()) == false){
				//System.out.println("ID: " + Integer.toString(n.getId()));
				//System.out.println("Distance: " + Integer.toString(n.getLengthToNode()));
				
				int nId = n.getId();
				int lengthToN = n.getLengthToNode();
				paths.put(n.getId(), lengthToN) ;
				if (this.incomingEdges.get(nId) != null){
					for (int i = 0; i < this.incomingEdges.get(nId).size(); i++) {
						int outEdge = this.incomingEdges.get(nId).get(i);
						int nIDLoc = Arrays.binarySearch( this.graph[outEdge], nId);
						int edgeLength = this.edges[outEdge][nIDLoc];
						int l = lengthToN + edgeLength;
						DijkstraNode d;
						if ( accessIntoHeap.containsKey(outEdge) ){
							d = accessIntoHeap.get(outEdge);
							if ( d.getLengthToNode() < l ){
								d.setLengthToNode(l);
							}
						}
						else{
							d = new DijkstraNode(outEdge, l);
						}
						pq.add(d);
	
					}
					//System.out.println();
				}
			}
		}

		int[] p = new int[this.graph.length];
		for (int a: paths.keySet()) {
			p[a] = paths.get(a); 
		}
		System.out.println( Arrays.toString(p) );
		return p;
	}
	
	public Matrix getAlphaFromSampledData(int[][] samples, Matrix[] alphaKStates){
		Matrix m = new Matrix( intArrayToDouble(samples) );
		Matrix durations = new Matrix( new double[][]{m.getArrayCopy()[0]} );
		Matrix distances = new Matrix( new double[][]{m.getArrayCopy()[1]} );
		
		double[][] a = new double[samples[0].length][alphaKStates[0].getArray().length];
		for (int i = 0; i < samples[0].length; i++) {
			int k = (int) durations.get(0, i);
			double[] aK = alphaKStates[k].getArrayCopy()[0];
			a[i] = aK;
		}
		Matrix A = new Matrix(a);
		Matrix AT = A.transpose();
		Matrix inverseForRegression = (AT.times(A)).inverse();
		
		Matrix theta = inverseForRegression.times( AT.times(distances.transpose()) ); //x = (At*A)^-1*AtD
		return A.times(theta);
	}
	
	public void performanceDistanceErrorComputations(Matrix Atheta, double[][] trueDistanceKAhead, int[][] samples){
		Matrix p = new Matrix( new double[][]{prior} ).transpose();
		Matrix td = new Matrix(trueDistanceKAhead);
		Matrix trueAverageDistance = td.times( p );	//Function of times
		
		Matrix s = new Matrix( intArrayToDouble(samples) );
		Matrix n = new Matrix( new double[][]{ {1,0} } );
		s = n.times(s);
		double[] distances = new double[s.getArray()[0].length];
		for (int i = 0; i < s.getArray()[0].length; i++) {
			distances[i] = trueAverageDistance.get( (int) s.getArray()[0][i], 0);
		}
		Matrix d = new Matrix(new double[][]{ distances }).transpose();
		
		Matrix error = Atheta.minus( d );

		d.transpose().print(5, 5);
		Atheta.transpose().print(5, 5);
		System.out.println( error.norm1() / error.getArray().length);
		//System.out.println( d.plus(Atheta).norm1() / (error.getArray().length*2) );
	}
	
	public int[][] createObservationDistanceSamples(int[] shortestPaths, int maxObservation, int samples){
		int[][] s = new int[2][samples];
		Random random = new Random();
		for (int i = 0; i < samples; i++) {
			int r = random.nextInt(maxObservation); 
			int d = sampleDistance(shortestPaths, r, random , prior);
			s[0][i] = r;
			s[1][i] = d;
		}
		
		return s;
	}

	private int sampleDistance(int[] shortestPaths, int r, Random random, double[] initialStates) {
		while(true){
			int state = sampleState(initialStates, random);
			//System.out.println("Initial State is: " + Integer.toString(state));
			int d = r;
			while(true){
				int desiredNextStateIndex = sampleState(transitions[state], random);
				if (desiredNextStateIndex == 0 && graph[state][0] == 0){
					break;
				}
				else if( edges[state][desiredNextStateIndex] > d){
					return shortestPaths[state] + d;
				}
				else{
					d = d - edges[state][desiredNextStateIndex]; 
					state = graph[state][desiredNextStateIndex];
				}
			}
		}
	}
	
	public int sampleState(double[] d, Random random){
		double[] c = super.cumulativeSum(d);
		double r = random.nextDouble();
		
		int index = Arrays.binarySearch(c, r);
		if (index < 0 ){
			index = index*-1 - 1;
		}

		return index;
	}
	
	public double[][] dynamicallyDetermineTrueDistanceKAhead(int[] shortestPaths, int maxK){
		double[][] distanceStorage = new double[maxK][this.graph.length];		
		
		for (int i = 0; i < maxK; i++) {
			for (int j = 1; j < this.graph.length; j++) {
				double[] tNoDoors = convertTransitionsToNonDoorTransitions(j);
				if (i == 0){
					distanceStorage[i][j] = shortestPaths[j];
				}
				else{
					for (int j2 = 0; j2 < graph[j].length; j2++) {
						int n = graph[j][j2];
						int d = i - edges[j][j2];
						if(d >= 0){
							distanceStorage[i][j] += distanceStorage[d][n]*tNoDoors[j2];
						}
						else{
							distanceStorage[i][j] += (distanceStorage[0][n] + -1*d)*tNoDoors[j2];
						}
					}
				}
				
			}
		}
		//System.out.println( Arrays.toString(distanceStorage[300]) );
		return distanceStorage;
	}
	
	public double[] convertTransitionsToNonDoorTransitions(int fromIndex){
		if(this.graph[fromIndex][0] == 0){
			double l = transitions[fromIndex].length;
			double[] t = new double[ (int) l];
			for (int i = 1; i < l; i++) {
				t[i] = transitions[fromIndex][i] * (l / (l-1.0));
			}
			return t;
		}
		else{
			return transitions[fromIndex];
		}
	}
	
	public static double[][] intArrayToDouble(int[][] I){
		double[][] r = new double[I.length][I[0].length];
		for (int j = 0; j < I.length; j++) {
			for (int j2 = 0; j2 < I[j].length; j2++) {
				//System.out.println( Arrays.toString(I[j]) );
				r[j][j2] = (double) I[j][j2];
			}
		}
		return r;
	}

}
