package hmm_sim;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Random;

import Jama.Matrix;

public class rawHMM extends Environment{
	
	private Matrix T;
	private Matrix O;
	private Matrix P;
	private Matrix E;
	public static Random random = new Random();
	private int automatonStates;
	
	public static void main(String[] args){
	
	}
	
	public rawHMM(String workingFolder, String description, int desiredHankelSize, Matrix T, Matrix O, Matrix P, Matrix E){
		super(workingFolder, description, desiredHankelSize);
		
		this.T = T;
		this.O = O;
		this.P = P;
		this.E = E;
		this.automatonStates = P.getArray()[0].length;
	}
	
	public static rawHMM make2StateHMM(String workingFolder){
		double[][] p = { {1}, {0}};
		double[][] t = { {0.5,0.45}, {0.3,0.67} };
		double[][] o = { {0,1}, {0,1} };
		double[][] e = { {0.05}, {0.03} };
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		Matrix E = new Matrix( e );
		
		int hSize = 100;
		rawHMM h = new rawHMM(workingFolder, "2_State_HMM",hSize,T,O,P,E);	
		return h;
	}
	
	public static rawHMM makeLabyrinth(String workingFolder, int loop1, int loop2 , double selfTransitionP, int hSize, double exit1p, double exit2p){
		boolean debug = false;
		
		int states = loop1 + loop2 - 1;
		HashMap<Integer, Double> termStates = new HashMap<Integer, Double>();
		int door1 = 0;
		int door2 = loop1/2 + loop2/2;
		termStates.put(door1, exit1p);
		termStates.put(door2, exit2p);
		
		HashMap<Integer, int[]> changeTo = new HashMap<Integer, int[]>();
		changeTo.put(loop1/2, new int[]{loop1/2 + 1, loop2 + loop1/2} );
		changeTo.put(states-1, new int[]{0});
		changeTo.put(loop2 + loop1/2 - 1, new int[]{loop1/2});

		double[][] p = new double[states][1];
		p[0][0] = 1;
		double[][] t = new double[states][states];
		double[][] e = new double[states][1];
		int[] v;
		for (int i = 0; i < states; i++) {
			if (changeTo.containsKey(i) ){	
			 	v = changeTo.get(i);
			 	for(int c=0;c<v.length;c++){
				   t[v[c]][i] = 1.0/v.length;
			 	}
			} 
			else if(termStates.containsKey(i)){
				t[i+1][i] = 1 - termStates.get(i);
				e[i][0] = termStates.get(i);
			}
			else{
				t[i+1][i] = 1-selfTransitionP;
				t[i][i] = selfTransitionP;
			} 
		}

		double[][] o = new double[states][states];
		for (int i = 0; i < states; i++) {
			o[i][i] = 1;
		}
		
		if (debug){
			System.out.println("Door1: ");
			System.out.println(door1);
			System.out.println("Door2: ");
			System.out.println(door2);
			System.out.println("From");
			System.out.println(loop1/2);
			System.out.println("To");
			System.out.println(loop1/2+1);
			System.out.println(loop1/2 + loop2);
			System.out.println("End");
			System.out.println(states-1);
			System.out.println("From");
			System.out.println(loop2 + loop1/2 - 1);
			System.out.println("To");
			System.out.println(loop1/2);
				
		}
		Matrix P = new Matrix( p ).transpose();
		Matrix T = new Matrix( t ).transpose();
		Matrix O = new Matrix( o );
		Matrix E = new Matrix( e );
		
		//String d = Integer.toString(loop1) + "_" + Integer.toString(loop2) + "_Toy_Labyrinth";
		rawHMM l = new rawHMM(workingFolder, workingFolder, hSize, T, O, P, E );
		
		return l;
		
	}
	
	/*public double[] generateEmpericalProbabilities(int samples) {
		int[] counts = new int[2*this.getDesiredHankelSize()];
		int  sequenceLength=0;
		
		for (int i = 0; i < samples; i++) {
			sequenceLength = this.generateSequenceStreakCount();	
			if(sequenceLength < counts.length){
				counts[sequenceLength]++;
			}
		}
		
		double[] probabilities = new double[counts.length];
		
		for (int i = 0; i < counts.length; i++) {
			probabilities[i] = ((double) counts[i])/samples;
		}
		
		return probabilities;
	}*/

	@Override
	public double[] generateTrueProbabilities() {
		Matrix P_True, S_True;
		
		Matrix Asigma = this.O.times(this.T);
		
		double[][] p = new double[2*this.getDesiredHankelSize()][this.getAutomatonStates()];
		double[][] s = new double[2*this.getDesiredHankelSize()][this.getAutomatonStates()];
		
		Matrix runningProductPrefixes = this.P;
		Matrix runningProductSuffixes = this.E;

		
		for (int i = 0; i < s.length; i++){
			p[i] = runningProductPrefixes.getArrayCopy()[0]; 
			s[i] = runningProductSuffixes.transpose().getArrayCopy()[0];
			runningProductPrefixes = runningProductPrefixes.times(Asigma);
			runningProductSuffixes = Asigma.times(runningProductSuffixes);
		}
		
		P_True = new Matrix(p);
		S_True = new Matrix(s).transpose();
		
		Matrix h = P_True.times(S_True);	

		return h.getArrayCopy()[0];
	}
	
	public int generateSequenceStreakCount(){
		
		int hiddenState = rawHMM.generateState( P.getArrayCopy()[0] );		
		int c = 0;
		double[] statesPossibilities;
		int states = P.getArrayCopy()[0].length;
		
		while(true){			
			statesPossibilities = T.getArrayCopy()[hiddenState];
			hiddenState = rawHMM.generateState( statesPossibilities );
			if( hiddenState < states){
				c++;
			}
			else{
				break;
			}
			
		} 

		return c;
	}
	
	public int[] generateSequence(int duration){
		int[] oSeq = new int[duration];
		
		int hiddenState = rawHMM.generateState( P.getArrayCopy()[0] );		
		
		for(int t=0;t<duration;t++){
			hiddenState = rawHMM.generateState( T.getArrayCopy()[hiddenState] );
			oSeq[t] = rawHMM.generateState( O.getArrayCopy()[hiddenState] );
		} 
		
		return oSeq;
	}
	
	public static int generateState(double[] stateProbabilities){ 						
		int l = stateProbabilities.length; 
		double[] cumulativeSum = new double[l];
		
		cumulativeSum[0] = stateProbabilities[0];
		for (int i = 1; i<l;i++){
			cumulativeSum[i] = cumulativeSum[i-1] + stateProbabilities[i];
		}

		double r = random.nextDouble();
		
		int index = Arrays.binarySearch(cumulativeSum, r);

		if (index >= 0){
			return index;
		}
		else{
			return -1*(index + 1);
		}
		
	}
	
	public int getAutomatonStates(){
		return this.automatonStates;
	}
}
