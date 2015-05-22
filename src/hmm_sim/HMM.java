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
		double[][] t = { {0.5,0.45}, {0.3,0.65} };
		double[][] o = { {0,1}, {0,1} };
	
		double[][] e = { {0.05}, {0.05} };
		
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		
		Matrix E = new Matrix( e );
		
		HMM h = new HMM(T,O,P,E,2);	
		
		//testEMGaussian();
		//testHankel();
		//h.testBaumWelch();
		//h.testFullBaumWelch();
		
		HashMap<String, Matrix> emp = h.singledataSpectralEmperical(1000,20);
		HashMap<String, Matrix> tru = h.singledataSpectralTrue(10);
		//System.out.println( emp.toString() );
		Matrix H = tru.get("P").times(tru.get("S"));
		Matrix Hbar = emp.get("H");
		
		System.out.println( H.minus(Hbar).normF() ); 
		
		System.out.println(emp.toString());
		
		Matrix temp1, temp2, temp3, r;
		for (int i = 0; i < 5; i++) {
			temp1 = emp.get( Integer.toString( (int) Math.pow(2, i)) );
			temp2 = HelperFunctions.matrixPower( temp1 , 2);
			temp3 = emp.get( Integer.toString( (int) Math.pow(2, i+1) ));
			r = temp2.minus( temp3 ) ;	
			//temp2.print(5, 5);
			//temp3.print(5, 5);
			
			System.out.println(i);
			System.out.println(r.normF());
		}
		
	}
	
	public HashMap<String, Matrix> singledataSpectralEmperical(int samples, int maxsize){
		int[] counts = new int[maxsize];
		int  sequenceLength=0, j=0, c=0, total=0;
		
		for (int i = 0; i < samples; i++) {
			sequenceLength = generateSequenceStreakCount();	
			if(sequenceLength < maxsize){
				counts[sequenceLength]++;
				total++;
			}
			
		}
		double[] probabilities = new double[total];
		
		for (int i = 0; i < maxsize; i++) {
			j = maxsize - i - 1;
			//counts[j] += c;
			
			probabilities[j] = ((double) counts[j])/samples;
			
			//c+=counts[j];							
			// The variable c makes counting occurrences linear instead of quadratic
		}
		
		//System.out.println( Arrays.toString(probabilities) );
		
		return singleObservationHankel(probabilities, 10, 2, states );
	}
	
	public HashMap<String, Matrix> singledataSpectralTrue(int size){
		Matrix P_True, S_True;
		
		Matrix Asigma = O.times(T);
		
		double[][] p = new double[size][states];
		double[][] s = new double[size][states];
		
		Matrix runningProductPrefixes = P.transpose();
		Matrix runningProductSuffixes = E;

		for (int i = 0; i < size; i++){
			p[i] = runningProductPrefixes.getArrayCopy()[0]; 
			s[i] = runningProductSuffixes.transpose().getArrayCopy()[0];
			runningProductPrefixes = runningProductPrefixes.times(Asigma);
			runningProductSuffixes = Asigma.times(runningProductSuffixes);
		}
		
		P_True = new Matrix(p);
		S_True = new Matrix(s).transpose();
		//P_True.print(5,5);
		//S_True.print(5,5);		

		//Test to verify above
		/*
		double[] sigmas = new double[size];
		Matrix running = P.transpose();
		for (int i = 0; i < size; i++) {
			sigmas[i] = running.times(E).trace();
			running = running.times(Asigma);
		}
		
		System.out.println( Arrays.toString(sigmas) );
		System.out.println( Arrays.toString(P_True.times(S_True).getArrayCopy()[0]) );
		*/
		
		HashMap<String, Matrix>  data= new HashMap<String, Matrix>();
		data.put("P", P_True);
		data.put("S", S_True);
		
		return data;
	}

	public void testBaumWelch(){
		int seq[] = new int[5];
		for (int i = 0; i < seq.length; i++) {
			seq[i] = generateSequenceStreakCount();

		}
		baumWelch(2, 10 ,seq );
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
	
	private int generateSequenceStreakCount(){
		
		int hiddenState = HelperFunctions.generateState( P.getArrayCopy()[0] );		
		int c = 0;
		while(true){
			//System.out.println(hiddenState);
			hiddenState = HelperFunctions.generateState( T.getArrayCopy()[hiddenState] );
			if( hiddenState <= 1){
				c++;
			}
			else{
				break;
			}
		} 
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
		
	public static void baumWelch(int numStates, int numIterations, int durations[] ){
		Matrix P,T,E, Pdiag;
		
		double[] prior =  HelperFunctions.randomNormVector(numStates);
		double[][] transition = new double[numStates][numStates]; 
		double[] ending = new double[numStates];
		
		for (int s = 0; s < numStates; s++) {
			double[] v = HelperFunctions.randomNormVector(numStates + 1);
			transition[s] = Arrays.copyOfRange(v, 0, v.length-1);
			ending[s] = v[v.length-1];
		}
		
		P = new Matrix(prior,1);
		T = new Matrix(transition);
		E = new Matrix(ending,1).transpose();
		
		double[][] pd = new double[numStates][numStates];
		for (int i = 0; i < pd.length; i++) {
			pd[i][i] = prior[i];
		}
		
		Pdiag = new Matrix(pd);
		
		//Compute probability of being in state i given a duration sequence
		double[][] forward;
		Matrix m,r;
		for (int d: durations){
			 forward = new double[d][numStates];
			 m = HelperFunctions.matrixPower(T,d);	//Inefficient fix later
			 r = Pdiag.times(m).times(E);
			 System.out.println(d);
			 r.print(5,10);	 
		}

	}
	
	public static HashMap<String, Matrix> singleObservationHankel(double[] counts, int basisSize , int base, int numHiddenStates){
		
		//	Matrix full = buildHankel(counts, 0, 20);
		//	full.print(5, 5);
			
		Matrix H = buildHankel(counts, 0, basisSize);

		//H.print(5, 5);
		//H = truncateSVD(H, numHiddenStates);
		//H.print(5, 5);
		
		SingularValueDecomposition SVD = H.svd();
		Matrix pinv = SVD.getU().times(SVD.getS()).inverse();
		Matrix sinv = (SVD.getV().transpose()).inverse();
						
		ArrayList<Matrix> H_Matrices  = new ArrayList<Matrix>();
		HashMap<String, Matrix> A_Matrices  = new HashMap<String, Matrix>();
		
		int maxDigit = (int) Math.floor( Math.log(counts.length )/Math.log(base) ) - 1; 
		int freq;
		Matrix h;
		for (int l = 0; l < maxDigit; l++) {
			freq = (int) Math.pow(base,l);
			h = buildHankel(counts, freq, freq+basisSize);
			//H_Matrices.add( truncateSVD(h, numHiddenStates) );
			H_Matrices.add(h);
		}
		
		Matrix m;
		for (int i = 0; i < H_Matrices.size(); i++) {
			m = pinv.times(H_Matrices.get(i)).times( sinv );
			A_Matrices.put(Integer.toString( (int) Math.pow(2, i) ), m );
		}
		
		double[][] h_L = new double[basisSize][1];
		for (int i = 0; i < basisSize; i++) {
			h_L[i][0] = (double) counts[i];
		}
		Matrix h_LS = new Matrix( h_L ).transpose();
		Matrix h_PL = h_LS.transpose();
				
		Matrix alpha_0 = h_LS.times(sinv);
		Matrix alpha_inf = pinv.times(h_PL);
		
		//SVD.getU().times(SVD.getS()).times(alpha_inf).print(5, 5); //Tests that you get back h_L
		
		/*alpha_0.times(A_Matrices.get("1")).times(alpha_inf).print(5, 5);
		alpha_0.times(A_Matrices.get("2")).times(alpha_inf).print(5, 5);
		alpha_0.times(A_Matrices.get("3")).times(alpha_inf).print(5, 5);
		*/
		A_Matrices.put("H", H);

		return A_Matrices;
	}
	
	public static Matrix truncateSVD(Matrix H, int nStates){
	
		SingularValueDecomposition svd = H.svd();
	    Matrix U = svd.getU();
	    Matrix S = svd.getS();
	    Matrix V = svd.getV();
	    
	    double[][] utemp = U.getArrayCopy();
	    double[][] utrunc = new double[utemp.length][nStates];
	    for (int i = 0; i < utrunc.length; i++) {
			utrunc[i] = utemp[i];
		}
	    
	    Matrix Utrunc = new Matrix(utrunc);
	    
	    double[][] stemp = S.getArrayCopy();
	    double[][] strunc = new double[stemp.length][nStates];
	    for (int i = 0; i < strunc.length; i++) {
			strunc[i] = stemp[i];
		}
	    
	    Matrix Strunc = new Matrix(strunc);
	    
	    double[][] vtemp = V.transpose().getArrayCopy();
	    double[][] vtrunc = new double[utemp.length][nStates];
	    for (int i = 0; i < vtrunc.length; i++) {
			vtrunc[i] = vtemp[i];
		}
	    
	    Matrix Vtrunc = new Matrix(vtrunc).transpose();
	    
	    return Utrunc.times(Strunc).times(Vtrunc);
	    
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
	
	
	//Not working and development halted because full HMM not needed for project
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
	
  }

	

