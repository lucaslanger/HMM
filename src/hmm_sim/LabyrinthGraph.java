package hmm_sim;

import java.awt.AlphaComposite;
import java.awt.image.SampleModel;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.PriorityQueue;
import java.util.Random;

import Jama.Matrix;
import Jama.SingularValueDecomposition;


public class LabyrinthGraph extends Environment{
	
	private int[][] graph;
	private int[][] edges;
	private double[][] transitions; 
	private double[] prior;
	private int stretchFactor;
	private int key;
	public boolean verbose;
	private HashMap<Integer, ArrayList<Integer>> incomingEdges;
	private HashMap<Integer, Integer> loopLengths;
	public HashMap<Integer, Integer> getLoopLengths() {
		return loopLengths;
	}

	public SequenceOfSymbols[] getData() {
		return data;
	}

	private SequenceOfSymbols[] data;
	
	
	public LabyrinthGraph(String workingFolder, int desiredHankelSize, int[][] graph, int[][] edges, double[][] transitions, double[] prior, int stretchFactor, int key, boolean verbose){
		super(workingFolder, workingFolder, desiredHankelSize);
		this.graph = graph;
		this.edges = edges;
		this.transitions = transitions;
		this.prior = prior;
		this.key = key;
		this.verbose = verbose;
		if(verbose){
			System.out.println("Key location is: ");
			System.out.println(key);
		}
		this.stretchFactor = stretchFactor;
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
		
		if (verbose){
			System.out.println("Probabilities:");
			r.transpose().print(5, 5);
			System.out.println("Sum");
			System.out.println(LabyrinthGraph.sumArray(r.transpose().getArrayCopy()[0]));
			System.out.println();
		}
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
	
	public static LabyrinthGraph pacMan(String workingFolder, int desiredHankelSize, int stretchFactor, int key, boolean verbose){
			
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

		 LabyrinthGraph l = new LabyrinthGraph(workingFolder, desiredHankelSize, graph, edges, transitions, prior, stretchFactor,key, verbose);
		 return l;
	}
	
	public static LabyrinthGraph multipleObservationDoubleLoop(String workingFolder, int desiredHankelSize, int numberOfTrajectories){
		int[][] graph = new int[][]{
				{1,2},
				{0},
				{0}
		};
		
		int[][] edges = new int[][]{
				{1,2},
				{1},
				{2}
		};
		
		double[][] transitions = new double[][]{
				{0.5,0.5},
				{1},
				{1}
		};
		
		int[][] wallColors = new int[][]{
				{1,2},
				{1},
				{2}
		};
		
		double[] prior = new double[]{1,0,0};
		
		LabyrinthGraph l = new LabyrinthGraph(workingFolder, desiredHankelSize, graph, edges, transitions, prior, 1, 1, false);
		HashMap<Integer, Integer> loopLengths = new HashMap<Integer, Integer>();
		int firstWallLength = wallColors[0][0] + wallColors[1][0];
		int secondWallLength = wallColors[0][1] + wallColors[2][0];
		loopLengths.put(1, firstWallLength);
		loopLengths.put(2, secondWallLength);
		
		l.loopLengths = loopLengths;
		
		l.data = l.generateSequencesMO(wallColors, numberOfTrajectories, desiredHankelSize, 0);
		return l;
	}
	
	// NOT GENERAL AT ALL, QUICK HACK TO OBTAIN TRUE PROBABILITY FOR SIMPLE ENVIRONMENT
	public double determineRealProbabilityOfSequenceDoubleLoop(SequenceOfSymbols seq){
		SequenceOfSymbols s = seq;
		double probability = 1;
		while(s.getSequence().equals("") == false){
			 SequenceOfSymbols fs = s.getFirstStreak();
			 String sym = fs.getSymbolFromString();
			 int duration = fs.getStreakFromString();
			 
			 int lS = this.loopLengths.get( Integer.parseInt(sym) );
			 int dmodL = duration % lS;
			 int multiples = (int) Math.floor(duration/lS);
			 
			 probability = probability * Math.pow(0.5, multiples);
			 
			 if( dmodL == 0){
				 if (fs.rawStringLength() == s.rawStringLength() ){
					 return probability;
				 }
				 else{
					 s = s.substring(fs.rawStringLength()+1, s.rawStringLength());
				 }
			} 
			 
			 else if(dmodL != 0 && fs.rawStringLength() == s.rawStringLength() ){
				 return 0.5*probability;
			 }

			 else{
				 return 0;
			 }
			 
		}
		return probability;
		
	}
	
	public SequenceOfSymbols[] generateSequencesMO(int[][] wallColors, int numberOfTrajectories, int maxLength, int startingLocation){
		SequenceOfSymbols[] seqs = new SequenceOfSymbols[numberOfTrajectories];
		for (int i = 0; i < numberOfTrajectories; i++) {
			seqs[i] = generateMultipleObservationSequence( wallColors, startingLocation, maxLength);
		}
		return seqs;
	}
	
	private SequenceOfSymbols generateMultipleObservationSequence(int[][] wallColors, int startingLocation , int maxLength) {
		SequenceOfSymbols s = new SequenceOfSymbols("");
		Random r = new Random();
		int trajLength = r.nextInt(maxLength);
		int lengthTravelled = 0;
		
		int currentState = startingLocation;
		while(lengthTravelled < trajLength){
			double[] nextStates = this.transitions[currentState];
			int nextStateIndex = sampleState( nextStates, r.nextDouble());
			int nextState = graph[currentState][nextStateIndex];
			int travelled = this.edges[currentState][nextStateIndex];
			int wallColor = wallColors[currentState][nextStateIndex];
			
			String a = Integer.toString(wallColor) + ":" + Integer.toString(travelled);
			SequenceOfSymbols seq = new SequenceOfSymbols(a);
			
			s = SequenceOfSymbols.concatenateSymbols(s, seq);
			
			lengthTravelled += travelled;
			currentState = nextState;
		}
		return s;
		
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

	public static LabyrinthGraph testLabyrinth(String workingFolder, int desiredHankelSize, int stretchFactor, boolean verbose){
		
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
		 
		 int key = 2;

		 LabyrinthGraph l = new LabyrinthGraph(workingFolder, desiredHankelSize, graph, edges, transitions, prior, stretchFactor, key, verbose);
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
		LabyrinthGraph.testLabyrinth(wf,  d, 1, false);
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
		//System.out.println( Arrays.toString(p) );
		return p;
	}
	
	public Matrix getAlphaFromSampledData(HashMap<String, int[]> trainingSamples, Matrix[] alphaKStates){

		Matrix distances = new Matrix( new double[][]{ intArrayToDouble(trainingSamples.get("Distances") ) } );
		Matrix durations = new Matrix( new double[][]{ intArrayToDouble(trainingSamples.get("Durations")) } );
			
		int numSamples = trainingSamples.get("Distances").length;
		int modelSize = alphaKStates[0].getArrayCopy().length;
		
		double[][] a = new double[numSamples][modelSize];
		for (int i = 0; i < numSamples; i++) {
			int k = (int) durations.get(0, i);
			double[] aK = alphaKStates[k].getArrayCopy()[0];
			a[i] = aK;
		}
		
		Matrix A = new Matrix(a);
		Matrix AT = A.transpose();
		Matrix ATA = AT.times(A);
				
		Matrix ATAinverse = svdInverse(ATA, false);
		//testInverse();

		Matrix theta = ATAinverse.times( AT.times(distances.transpose()) ); //x = (At*A)^-1*AtD
		
		return theta;
	}
	
	public double norm2Custom(double[] d){
		double r = 0;
		for (int i = 0; i < d.length; i++) {
			r += Math.pow(d[i], 2);
		}
		return Math.sqrt(r);
	}
	
	public double norm1Custom(double[] d){
		double r = 0;
		for (int i = 0; i < d.length; i++) {
			r += Math.abs(d[i]);
		}
		return r;
	}
	
	public void testInverse(){
		double[][] d = { { 1,3,2}, {5,3,1}, {2,0,0}};
		Matrix D = new Matrix(d);
		D.times( D.inverse() ).print(5, 5);
		D.times( svdInverse(D, false) ).print(5, 5);
	}
	
	public Matrix svdInverse(Matrix m, boolean debug){
		SingularValueDecomposition svd = m.svd();
		if (debug == true){
			//svd.getS().print(5, 5);
			//svd.getU().print(5, 5);
			//svd.getV().print(5, 5);
			svd.getV().transpose().times( svd.getV()).print(5, 5);
			svd.getU().transpose().times( svd.getU()).print(5, 5);
			pseudoInvDiagonalKillLowSingularValues(svd.getS()).print(5, 5);
		}
		return svd.getV().times( pseudoInvDiagonalKillLowSingularValues(svd.getS() ) ).times(svd.getU().transpose() );
	}
	
	public Matrix pseudoInvDiagonalKillLowSingularValues(Matrix m){
		double[][] a = m.getArrayCopy();
		for (int i = 0; i < a.length; i++) {
			if (a[i][i] > 0.0001){
				a[i][i] = 1/a[i][i];
			}
			else{
				a[i][i] = 0;
			}
		}
		return new Matrix(a);
	}
	/*
	public double performanceDistanceErrorComputations(Matrix Atheta, double[][] trueDistanceKAhead, int[][] samples){
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
		
		return error.norm1()/ error.getArray().length;
	}*/
	
	public HashMap<String, int[]> createObservationDistanceSamples(int[] shortestPaths, int maxObservation, int samples){
		HashMap<String, int[]> s = new HashMap<String, int[]>();
		Random random = new Random();
		int[] durations = new int[samples];
		int[] distances = new int[samples];
		for (int i = 0; i < samples; i++) {
			int r = random.nextInt(maxObservation); 
			int d = sampleDistance(shortestPaths, r, random , prior);
			durations[i] = r;
			distances[i] = d;
		}
		s.put("Distances", distances);
		s.put("Durations", durations);
		
		return s;
	}

	private int sampleDistance(int[] shortestPaths, int r, Random random, double[] initialStates) {
		while(true){
			int state = sampleState(initialStates, random.nextDouble());
			//System.out.println("Initial State is: " + Integer.toString(state));
			int d = r;
			while(true){
				int desiredNextStateIndex = sampleState(transitions[state], random.nextDouble());
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
	
	public static int sampleState(double[] a, double r){
		double[] c = cumulativeSum(a);
		
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
	
	public static double[] intArrayToDouble(int[] I){
		double[] r = new double[I.length];
		for (int j = 0; j < I.length; j++) {
			r[j] = (double) I[j];
		}
		return r;
	}
	
	/*
	public Matrix[] getTrueAlphaKs(int k, int maxK){		
		int trueStateSize = graph.length - 1;
		int indexToExpandOn = graph.length; 
		ArrayList<Integer> indexMarkers = new ArrayList<Integer>();
		HashMap<IntegerPair, IntegerPair> pairsSeen = new HashMap<IntegerPair, IntegerPair>();
		
		for (int i = 0; i < edges.length; i++) {
			for (int j = 0; j < edges[i].length; j++) {
				IntegerPair temp = new IntegerPair(i, graph[i][j]);
				if (i ==0 || graph[i][j] == 0 || (pairsSeen.containsKey(temp) ) ){
					continue;
				}
				else{
					trueStateSize += edges[i][j] - 1;
					IntegerPair pair = new IntegerPair(i, graph[i][j]);
					IntegerPair range = new IntegerPair(indexToExpandOn, indexToExpandOn + edges[i][j] - 1);
					pairsSeen.put(pair, range);
					indexToExpandOn += edges[i][j] - 1;
				}
			}
		}
		
		int[] distribution = new int[trueStateSize];
		for (int i = 0; i < prior.length; i++) {
			off by 1 because of 0 state?
			distribution[i] = prior[i];
			
			FIX LABYRINTH TO BE A TRUE HMM
		}
		
		Matrix[] maxKs = new Matrix[maxK];
		for (int i = 1; i < maxK; i++) {
			int[] nextDistribution = new int[trueStateSize];
			for (int j = 0; j < trueStateSize; j++) {
				if (j < graph.length){
					double[] a = convertTransitionsToNonDoorTransitions(j);
					for (int out = 0; out < t.length; out++) {
						 
					} 
				}
				else{
					nextDistribution[i] = nextDistribution[i-1];
				}
			}
		}
		
		return maxKs;
	}
	*/
}
