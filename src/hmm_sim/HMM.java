//HMM model written by Lucas Langer 

package hmm_sim;

import Jama.*;

import java.util.*;

import hmm_sim.HelperFunctions;

public class HMM {
	
	// T --> Transition, O --> Observation, P --> Prior
	private int states;
	
	private Matrix T;
	private Matrix O;
	private Matrix P;
	private Matrix E;

	public static Random random = new Random();
	
	public static void main(String[] args){
		
		/* ONLY FOR TRUEHMM
		double[][] p = { {0}, {1}};
		double[][] t = { {0.55,0.45}, {0.35,0.65} };
		double[][] o = { {0,1}, {0,1} };
		*/
		
		double[][] p = { {0}, {1}};
		double[][] t = { {0.5,0.45}, {0.3,0.67} };
		double[][] o = { {0,1}, {0,1} };
	
		double[][] e = { {0.05}, {0.03} };
		
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		
		Matrix E = new Matrix( e );
		
		HMM h = new HMM(T,O,P,E,2);	
		
		//testEMGaussian();
		//testHankel();
		//h.testBaumWelch();
		//h.testFullBaumWelch();
		
		HashMap<String, Matrix> emp = h.singledataSpectralEmperical(100,100000,5);
		HashMap<String, Matrix> tru = h.singledataSpectralTrue(100,5);

		
		Matrix H = tru.get("H");
		Matrix Hbar = emp.get("H");	
		
		Matrix temp1, temp2, r;
		int pow;
		for (int i = 0; i < emp.get("max").norm1()-1; i++) {
			pow = (int) Math.pow(2, i);
			temp1 = emp.get( Integer.toString( pow ) );
			temp1 = HelperFunctions.matrixPower( temp1 , 2);
			temp2 = emp.get( Integer.toString(pow*2) );
			r = temp2.minus( temp1 ) ;	
			
		}
		
	}
		
	public HashMap<String, Matrix> singledataSpectralEmperical(int size, int samples, int basisSize){
		int[] counts = new int[2*size];
		int  sequenceLength=0, j=0;
		
		for (int i = 0; i < samples; i++) {
			sequenceLength = generateSequenceStreakCount();	
			if(sequenceLength < size){
				counts[sequenceLength]++;
			}
		}

		
		double[] probabilities = new double[counts.length];
		
		for (int i = 0; i < counts.length; i++) {
			probabilities[i] = ((double) counts[i])/samples;
		}
		
		return singleObservationHankel(probabilities, basisSize, 2, states );
	}
	
	public HashMap<String, Matrix> singledataSpectralTrue(int size, int basisSize){
		Matrix P_True, S_True;
		
		Matrix Asigma = O.times(T);
		
		double[][] p = new double[2*size][states];
		double[][] s = new double[2*size][states];
		
		//System.out.println(states);
		//Asigma.print(5, 5);
		
		Matrix runningProductPrefixes = P.transpose();
		Matrix runningProductSuffixes = E;

		for (int i = 0; i < s.length; i++){
			p[i] = runningProductPrefixes.getArrayCopy()[0]; 
			s[i] = runningProductSuffixes.transpose().getArrayCopy()[0];
			runningProductPrefixes = runningProductPrefixes.times(Asigma);
			runningProductSuffixes = Asigma.times(runningProductSuffixes);
		}
		
		P_True = new Matrix(p);
		S_True = new Matrix(s).transpose();
		
		Matrix h = P_True.times(S_True);	

		return singleObservationHankel(h.getArrayCopy()[0], basisSize, 2, states);
	}

	public HMM(Matrix T, Matrix O, Matrix P, Matrix E, int ns){
		this.T = T;
		this.O = O;
		this.P = P;
		this.E = E;
		this.states = ns;
	}
	
	public HMM(Matrix T, Matrix O, Matrix P,  int ns){
		this.T = T;
		this.O = O;
		this.P = P;
		this.states = ns;
	}
	
	public static HashMap<String, Matrix> singleObservationHankel(double[] counts, int basisSize , int base, int numHiddenStates){
			
		/*Test of concatenation method 
		Matrix cH = HelperFunctions.buildConcatHankel(counts, basisSize);
		HashMap<String, Matrix> SVDc = HelperFunctions.truncateSVD(cH, numHiddenStates);
		Matrix pinvc = (SVDc.get("U").times(SVDc.get("S"))).inverse();
		Matrix sinvc = SVDc.get("VT").transpose();
		System.out.println("Concat pinv");
		pinvc.print(5, 5);
		*/
		
		Matrix H = buildHankel(counts, 0, basisSize);

		HashMap<String, Matrix> SVD = HelperFunctions.truncateSVD(H, numHiddenStates);
		H = SVD.get("U").times(SVD.get("S")).times(SVD.get("VT"));	
		Matrix pinv = (SVD.get("U").times(SVD.get("S"))).inverse();
		Matrix sinv = (SVD.get("VT")).transpose();
						
		ArrayList<Matrix> H_Matrices  = new ArrayList<Matrix>();
		HashMap<String, Matrix> returnData  = new HashMap<String, Matrix>();
		
		int maxDigit = (int) Math.floor((Math.log( (counts.length/2) - basisSize)/Math.log(base))) ; //Too low fix later to allow higher powers
		int freq;
		Matrix h;
		for (int l = 0; l <= maxDigit; l++) {
			freq = (int) Math.pow(base,l);
			h = buildHankel(counts, freq, freq+basisSize);
			H_Matrices.add(h);
		}
		
		Matrix m;
		for (int i = 0; i < H_Matrices.size(); i++) {
			m = pinv.times(H_Matrices.get(i)).times( sinv );
			returnData.put(Integer.toString( (int) Math.pow(2, i) ), m );
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < basisSize; i++) {
			h_L[i][0] = (double) counts[i];
		}
		Matrix h_LS = new Matrix( h_L ).transpose();
		Matrix h_PL = h_LS.transpose();
				
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_PL);
				
		Matrix maX = new Matrix(new double[][]{{maxDigit}});
		
		returnData.put("H", H);
		returnData.put("max", maX);
		returnData.put("pinv", pinv);
		returnData.put("sinv", sinv);
		returnData.put("s_values", SVD.get("S"));
		returnData.put("U", SVD.get("U"));
		returnData.put("VT", SVD.get("VT"));
		returnData.put("a0", alpha_0);
		returnData.put("ainf", alpha_inf);

		return returnData;
	}

	public static Matrix buildHankel(double[] counts, int startingindex, int endingindex){
		int size = endingindex - startingindex;
		double[][] hankel = new double[size][size];
		
		int i,j;
		for (i = 0; i < size; i++) {
			for (j = 0; j < size; j++) {
				hankel[i][j] = counts[i+j+startingindex];
			}
		}
		
		return new Matrix(hankel);
	}
	
	//Development halted because full HMM not needed for project
	//Fix: Vectorize and debug non normal transition matrix
	
	public static void baumWelchFULLHMM(int numStates, int numIterations, int[] observations){
		double[][] prior =  { HelperFunctions.randomNormVector(numStates) };
		double[][] transition = new double[numStates][numStates]; 
		double[][] observ = new double[][]{ HelperFunctions.randomNormVector(2) , HelperFunctions.randomNormVector(2) };
		
		for (int i = 0; i < numStates; i++) {
			transition[i] = HelperFunctions.randomNormVector(numStates);
		}
			
		Matrix P,T,O,E;
		int duration = observations.length;
		double[][] forwardData = new double[ duration ][numStates];
		
		double[][][] transProbs = new double[duration-1][numStates][numStates];
		double[][] stateProbs = new double[duration][numStates];
		
		for (int c = 0; c < numIterations; c++) {
	
			P = new Matrix(prior);
			T = new Matrix(transition).transpose(); 
			O = new Matrix(observ);
			
			System.out.println("Run");
			P.print(numStates, numStates);
			T.print(numStates, numStates);
			O.print(numStates, numStates);
			
			int obs = 0;
			for (int i = 0; i < numStates; i++) {
				obs = observations[i];
				forwardData[0][i] = observ[obs][i] * prior[0][i];	
			}
			
			double tmp = 0;
			double tmp2 = 0;
			Matrix forw,tran;
			for (int i = 1; i < duration; i++) {
				for (int j = 0; j < numStates; j++) {
					tmp = 0;
					forw = new Matrix(forwardData[i-1],1).transpose();
					tran = new Matrix(transition[j], 1); 
					tmp = forw.times(tran).trace();
					
					tmp2 = 0;
					for (int p = 0; p < numStates; p++) {
						tmp += forwardData[i-1][p]*transition[j][p];	// from p to j
					}
					forwardData[i][j] = tmp*observ[obs][j];
				}
			}
			
			//System.out.println("Forward");
			//new Matrix(forwardData).print(5, 5);
			
			double[][] backwardData = new double[duration][numStates];
			for (int i = 0; i < numStates; i++) {
				backwardData[duration-1][i] = 1;
			}
			
			//Improvement: Turn the code below into matrix operations, numStates, observ states low so it doesn't matter
			int j;
			for (int i = 1; i < duration; i++) {
				for (int k = 0; k < numStates; k++) {
					tmp = 0;
					j = backwardData.length-1-i;
					for (int p = 0; p < numStates; p++) {
						if (observations[i] == 1){
							tmp += transition[p][k] * backwardData[j+1][p] * observ[0][p];
						}
						else{
							tmp += transition[p][k] * backwardData[j+1][p] * (1-observ[0][p]);
						}
					}
					if (tmp > 1) {
						//System.out.println(tmp);
					}
					
					backwardData[j][k] = tmp;
				}
			}
			
			//System.out.println("Backward");
			//new Matrix(backwardData).print(5, 5);
			
			double sum;
			for (int i = 0; i < duration; i++) {
				sum = 0;
				for (int k = 0; k < numStates; k++) {
					stateProbs[i][k] = forwardData[i][k]*backwardData[i][k];
					sum += stateProbs[i][k];
				}
				for (int k = 0; k < numStates; k++) {
					stateProbs[i][k] = stateProbs[i][k]/sum; 
				}
			}
			
			for (int i = 0; i < duration-1; i++) {
				sum = 0;
				for (int k = 0; k < numStates; k++) {
					for (int p = 0; p < numStates; p++) {
						if (observations[p] == 1){
							transProbs[i][k][p] = forwardData[i][k]*transition[p][k]*backwardData[i+1][p]*observ[0][p];
						}
						else{
							transProbs[i][k][p] = forwardData[i][k]*transition[p][k]*backwardData[i+1][p]*(1-observ[0][p]);
						}
						sum += transProbs[i][k][p];
					}
				}
				
				for (int k = 0; k < numStates; k++) {
					for (int p = 0; p < numStates; p++) {
						transProbs[i][k][p] = transProbs[i][k][p]/sum;
						if (transProbs[i][k][p] > 1) {
							System.out.println("Transprob");
							System.out.println(transProbs[i][k][p]);
						}
					}
				}
				
			}
			
			for (int i = 0; i < numStates; i++) {
				prior[0][i] = stateProbs[0][i];
			}
			
			double norm;
			for (int i = 0; i < numStates; i++) {
				for (int k = 0; k < numStates; k++) {
					norm = 0;
					for (int t = 0; t < duration-1; t++) {
						norm += stateProbs[t][i];
					}
					
					transition[i][k] = 0;
					for (int t= 0; t < duration-1; t++) {
						transition[i][k] += transProbs[t][i][k];
					}
					
					
					transition[i][k] = transition[i][k]/norm;
				}
			}
			
			double norm2, obscount;
			for (int i = 0; i < numStates; i++) {
				norm2 = 0;
				obscount = 0;
				for (int t = 0; t < duration; t++) {
					if( observations[t] == 1 ){
						obscount += stateProbs[t][i];
					}
					norm2 += stateProbs[t][i];
				}
				observ[0][i] = obscount/norm2;
			}
		}
		
		//Computation of P(Si, O1 .. Oi): Forward
		//Computation of P(Oi+1, ... On|Si) Backward
		//Now compute P(Si, O1...On) --> condition on Si, then you have a product
		//Transitions: P(Si,Si+1|O1...On), compute using info above and transitions guess
		//Special for Baum-Welch: since T's are fixed, compute sum over transitions from ij/ all possible transitions- average measure 
		//To compute O: taking sum t=1 to N #times Si=x*obs_frequency/ sum1toN obs_frequency
		//Idea above: Weighting observation likelyhood by likelyhood of Si=x there for t=1 to N
	}
	
	public static void testHankel(){
		Matrix Hbar = new Matrix( new double[][]{ {0,0.2,0.14}, {0.2,0.22,0.15}, {0.14,0.45,0.31} }).transpose();
		
		Matrix Ha = new Matrix(new double[][]{ {0.2,0.22,0.15},{0.22,0.19,0.13},{0.45,0.45,0.32} }).transpose();
		Matrix Hb = new Matrix(new double[][]{ {0.14,0.45,0.31}, {0.15,0.29,0.13}, {0.31,0.85,0.58} } ).transpose();
		Matrix hls = new Matrix(new double[][]{ {0, 0.2, 0.14} } );
		Matrix hpl = new Matrix(new double[][]{ {0, 0.2, 0.14} } ).transpose();
			
		Hbar.print(5,5);
		
		SingularValueDecomposition svd = Hbar.svd();
		Matrix p = svd.getU().times( svd.getS() );
		Matrix s = svd.getV().transpose();
		
		Matrix pinv = p.inverse();
		Matrix sinv = s.inverse();
		
		Matrix Aa = pinv.times(Ha).times(sinv); 
		Matrix Ab = pinv.times(Hb).times(sinv); 
		
		Matrix alpha0 = hls.times(sinv);	//alpha0 row
		Matrix alphainf = pinv.times(hpl);	//alphainf column
		
		
		Matrix test1 = alpha0.times(Aa).times(alphainf);
		Matrix test2 = alpha0.times(Ab).times(alphainf);
		Matrix test3 = alpha0.times(Aa).times(Ab).times(alphainf);
		Matrix test4 = alpha0.times(Ab).times(Aa).times(alphainf);
		Matrix test5 = alpha0.times(Aa).times(Aa).times(alphainf);
		Matrix test6 = alpha0.times(Ab).times(Ab).times(alphainf);
		
		test1.print(5,5);
		test2.print(5,5);
		test3.print(5,5);
		test4.print(5,5);
		test5.print(5,5);
		test6.print(5,5);
	}
		
	public void testFullBaumWelch(){
		int[] seq = generateSequence(20);
		//System.out.println(Arrays.toString(seq));
		baumWelchFULLHMM(2, 10, seq);
	}
	
	public static void testEMGaussian(){
		
		double c1 = random.nextDouble();
		double c2 = random.nextDouble();
		double c3 = random.nextDouble();
		
		System.out.println(c1);
		System.out.println(c2);
		System.out.println(c3);
		
		int N = 10;
		
		double[][] p = { {0}, {1}, {0} };
		double[][] t = { {0.3,0.3,0.4}, {0.3,0.3,0.4}, {0.3,0.3,0.4} }; 
		
		double[][] o = { 	HelperFunctions.generateBinomialVector(N,c1), 
							HelperFunctions.generateBinomialVector(N,c2), 
							HelperFunctions.generateBinomialVector(N,c3) 
						};
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		
		HMM h = new HMM(T,O,P,3);
		
		int[] int_sequence = h.generateSequence(100);
		double[] sequence = new double[100];
		for (int i = 0; i < sequence.length; i++) {
			sequence[i] = (double) int_sequence[i];
		}
				
		Set<String> s = HelperFunctions.expectationMaximizationGaussian(sequence, 3 , 50);
		
		System.out.println( s );
	}
	
	public int generateSequenceStreakCount(){
		
		int hiddenState = HelperFunctions.generateState( P.getArrayCopy()[0] );		
		int c = 0;
		double[] statesPossibilities;
		while(true){
			//System.out.println("HiddenState");
			//System.out.println(hiddenState);
			
			statesPossibilities = T.getArrayCopy()[hiddenState];
			//System.out.println(Arrays.toString(statesPossibilities));
			hiddenState = HelperFunctions.generateState( statesPossibilities );
			if( hiddenState < states){
				c++;
			}
			else{
				break;
			}
			
		} 
		//System.out.println("Streak");
		//System.out.println(c);
		return c;
	}
	
	public int[] generateSequence(int duration){
		int[] oSeq = new int[duration];
		
		int hiddenState = HelperFunctions.generateState( P.getArrayCopy()[0] );		
		
		for(int t=0;t<duration;t++){
			hiddenState = HelperFunctions.generateState( T.getArrayCopy()[hiddenState] );
			//System.out.println(hiddenState);
			oSeq[t] = HelperFunctions.generateState( O.getArrayCopy()[hiddenState] );
		} 
		
		return oSeq;
	}
	
  }

	

