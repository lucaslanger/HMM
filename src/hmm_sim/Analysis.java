package hmm_sim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import Jama.Matrix;

public class Analysis {
	
	public static final String pltFolder = "plotting/";
	
	private ArrayList<HashMap<String, Matrix>> empArray = new ArrayList<HashMap<String,Matrix>>();
	private HashMap<String, Matrix> tru;
	private HMM h;
	
	private int basisSize;
	private int hSize;
	
	public static void main(String[] args){
		Analysis a = new Analysis(100,30,1000,1000);
		
		a.compareSigmaError();
		System.out.println("Power Sigma Done");
		a.compareH_Hbar();
		System.out.println("H Hbar Done");
		a.compareASigmas();
		System.out.println("SigmaError Done");
		a.compareQueryErrors();
		System.out.println("Query Done");
		
	}
	
	public Analysis(int hSize, int basisSize, int trials, int trialSize){
		//this.h = makeHMM();
		this.h = makeLabyrinth();
		
		this.tru = h.singledataSpectralTrue(hSize, basisSize);
		
		this.tru.get("H").print(5,5); 
		
		this.hSize = hSize;
		this.basisSize = basisSize;
		
		for (int i = 0; i < trials; i++) {
			empArray.add(h.singledataSpectralEmperical(hSize, trialSize, basisSize));
		}
	}
	
	public void compareH_Hbar(){
		
		int[] sizes = new int[]{100,1000,10000};
		
		HashMap<String, Matrix> emp;
		
		Matrix H;
		Matrix Hbar;
		
		int trailsize = sizes.length;
		int repeats = 10;
		
		double[][] dataSize = new double[1][trailsize];
		double[][] error = new double[1][trailsize];
		double avgError, e;
	
		for (int i = 0; i < trailsize; i++) {
			avgError = 0;
			for (int j = 0; j < repeats; j++) {
				emp = h.singledataSpectralEmperical(hSize, sizes[i], basisSize);
				
				H = tru.get("H");
				Hbar = emp.get("H");
				
				e = H.minus(Hbar).normF();
				avgError += e;
			}
			avgError /= repeats;
			dataSize[0][i] = sizes[i];
			error[0][i] = avgError;
		}
			
		HelperFunctions.outputData(pltFolder + "True_Hankel_vs_Emperical", "X:#Data Seen Y:Fnorm","", dataSize, error );
		
	}
	
	public void compareASigmas(){
		HashMap<String, Matrix> emp;
		
		int maxQuery = (int) (tru.get("max").norm1()-1); 
		
		double[][] sigmaNumber = new double[1][maxQuery];
		double[][] errors = new double[1][maxQuery];
		
		Matrix h_sigma_true, h_sigma_exp, r;
		int pow;
		
		for (int j = 0; j < empArray.size(); j++) {
			emp = empArray.get(j);
			for (int i = 0; i < maxQuery; i++) {
				pow = (int) Math.pow(2, i);
				h_sigma_true = tru.get( Integer.toString(pow) );
				h_sigma_exp = emp.get( Integer.toString(pow) );
				r = h_sigma_true.minus( h_sigma_exp );
				errors[0][i] += r.normF();
			}
		}
		
		for (int i = 0; i < maxQuery; i++) {
			pow = (int) Math.pow(2, i);
			h_sigma_true = tru.get( Integer.toString(pow) );
			sigmaNumber[0][i] = pow;
			errors[0][i] /= (empArray.size()*h_sigma_true.normF());
		}
		
		HelperFunctions.outputData(pltFolder + "True_Ax_vs_Emperical_Ax", "X:Sigma Y:Relative Fnorm","", sigmaNumber, errors );
		// Add file containing error analysis for alphaInf and alpha0?
	}
	
	public void compareSigmaError(){
		int m = (int) (tru.get("max").norm1()-1); 
		
		double[][] sigmaNumber = new double[1][m];
		double[][] errors = new double[1][m];
		
		Matrix temp1, temp2,  r;
		int pow;
		
		for (int j = 0; j < empArray.size(); j++) {
				
			for (int i = 0; i < m; i++) {
				pow = (int) Math.pow(2, i);
				temp1 = empArray.get(j).get( Integer.toString( pow ) );
				temp1 = HelperFunctions.matrixPower( temp1 , 2);
				temp2 = empArray.get(j).get( Integer.toString(pow*2) );
				r = temp2.minus( temp1 ) ;	
				//r.print(5, 5);
			
				errors[0][i] += r.normF();
			}
			
		}
		
		Matrix h_sigma_true;
		for (int i = 0; i < m; i++) {
			pow = (int) Math.pow(2, i+1);
			h_sigma_true = tru.get( Integer.toString(pow) );
			sigmaNumber[0][i] = pow;
			errors[0][i] /= (empArray.size()*h_sigma_true.normF());
		}
		
		HelperFunctions.outputData(pltFolder + "(Ax)^2_v.s A(x^2)", "X:Sigma Y:Relative Fnorm","", sigmaNumber, errors );

	}
	
	public void compareQueryErrors(){
		int maxexp = (int) tru.get("max").norm1(); //same for tru by construction
		int maxquery = (int) Math.pow(2, maxexp);
		
		double[][] queries = new double[9][maxquery];
		double[][] errors = new double[9][maxquery];
		
		Matrix a0emp, ainfemp, empQF, empQB, empP, empProbQF, empProbQB, empProbP;
		Matrix truF, truB, a0tru, ainftru, truProbQF, truProbQB;
		
		HashMap<String, Matrix> emp;
		for (int i = 0; i < maxquery ; i++) {
			truF = HelperFunctions.matrixQuery(tru, i, 2, true);
			truB = HelperFunctions.matrixQuery(tru, i, 2, false);
			a0tru = tru.get("a0");
			ainftru = tru.get("ainf");
			truProbQF = a0tru.times(truF).times(ainftru);
			truProbQB = a0tru.times(truB).times(ainftru);
			//truProbF.minus(truProbB).print(5, 5); //Always 0 which makes sense
			
			/*
			System.out.println("Query, Prob");
			System.out.println(i);
			truProbQF.print(5,5);
			System.out.println(truProbQF.get(0,0));
			*/
			
			for (int j = 0; j < empArray.size(); j++) {	
				emp = empArray.get(j);
		
				empQF = HelperFunctions.matrixQuery(emp, i, 2, true);
				empQB = HelperFunctions.matrixQuery(emp, i, 2, false);

				empP = HelperFunctions.matrixPower(emp.get("1"), i);										//inefficient, if slow optimize later
				
				a0emp = emp.get("a0");
				ainfemp = emp.get("ainf");
				
				empProbQF = a0emp.times(empQF).times(ainfemp);
				empProbQB = a0emp.times(empQB).times(ainfemp);
				empProbP = a0emp.times(empP).times(ainfemp);
				
				/*
				System.out.println("Forward");
				empQF.print(5, 5);
				System.out.println("Backward");
				empQB.print(5, 5);
				System.out.println("a0");
				a0emp.print(5, 5);
				ainfemp.print(5, 5);
				empProbQF.minus(empProbQB).print(5, 5);
				*/
				
				
				errors[0][i] += Math.abs(truProbQF.minus(empProbQF).get(0,0));
				errors[1][i] += truProbQF.minus(empProbQF).get(0,0);

				errors[2][i] += Math.abs(truProbQF.minus(empProbP).get(0,0));
				errors[3][i] += truProbQF.minus(empProbP).get(0,0);
								
				errors[4][i] += Math.abs(empProbQF.minus(empProbQB).get(0,0));
				errors[5][i] += empProbQF.minus(empProbQB).get(0,0);
				
				errors[6][i] += empQF.minus(empQB).normF();
				
				errors[7][i] += Math.abs(truProbQF.minus(empProbQF).get(0,0))/ Math.max( Math.abs(truProbQF.get(0,0)), Math.abs(empProbQF.get(0,0)) );
				errors[8][i] += Math.abs(truProbQF.minus(empProbP).get(0,0))/ Math.max( Math.abs(truProbQF.get(0,0)), Math.abs(empProbP.get(0,0)) );
			}
			for (int c = 0; c < 9; c++) {
				errors[c][i] /= (empArray.size());
				queries[c][i] = i;
			}

			
		}
		
		HelperFunctions.outputData(pltFolder + "Query_Errors_Base", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,0,2), Arrays.copyOfRange(errors,0,2) );
		HelperFunctions.outputData(pltFolder + "Query_Errors_Naive", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,2,4), Arrays.copyOfRange(errors,2,4) );
		HelperFunctions.outputData(pltFolder + "Non-Comm_Query_Error", "X:Sigma Y:a0(A16A1-A1A16)aI","", Arrays.copyOfRange(queries,4,6), Arrays.copyOfRange(errors,4,6) );
		HelperFunctions.outputData(pltFolder + "Non-Comm_Matrix_Error", "X:Sigma Y:(A16A1-A1A16).Fnorm","", Arrays.copyOfRange(queries,6,7), Arrays.copyOfRange(errors,6,7) );
	
		double[][] ebase = Arrays.copyOfRange(errors,0, 1);
		double[][] enaive = Arrays.copyOfRange(errors, 2, 3);
		double[][] ejoint = new double[][]{ebase[0], enaive[0]};
		
		double[][] qbase = Arrays.copyOfRange(queries, 2, 3);
		double[][] qnaive = Arrays.copyOfRange(queries, 2, 3);
		double[][] qjoint = new double[][]{qbase[0], qnaive[0]};
		HelperFunctions.outputData(pltFolder + "QError_Base_vs_Naive", "X:Sigma Y:Absolute error","",qjoint,ejoint  );
		
		HelperFunctions.outputData(pltFolder + "QError_Rel_Base_vs_Naive", "X:Sigma Y:Absolute error","", qjoint, Arrays.copyOfRange(errors,7,9) );

	}
	
	public HMM makeHMM(){
		double[][] p = { {0.5}, {0.5}};
		double[][] t = { {0.5,0.45}, {0.3,0.67} };
		double[][] o = { {0,1}, {0,1} };
		double[][] e = { {0.05}, {0.03} };
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		Matrix E = new Matrix( e );
		
		HMM h = new HMM(T,O,P,E,2);	
		return h;
	}
	
	public HMM makeLabyrinth(){
		int states = 27;
		double selfTransitionP = 0.05;
		
		HashMap<Integer, Double> termStates = new HashMap<Integer, Double>();
		termStates.put(0, .5);
		termStates.put(13, .5);
		
		HashMap<Integer, int[]> changeTo = new HashMap<Integer, int[]>();
		changeTo.put(7, new int[]{8,20} );
		changeTo.put(26, new int[]{0});
		changeTo.put(19, new int[]{7});

		
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
		
		Matrix T = new Matrix( t ).transpose();
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		Matrix E = new Matrix( e );
		
		/*T.print(5, 5);
		E.print(5, 5);
		O.print(5, 5);
		P.print(5, 5);*/
		HMM l = new HMM(T, O, P, E, states);
		
		return l;
		
	}
	
	public void debugHComparisons(HashMap<String, Matrix> emp ){
		System.out.println("H error");
		emp.get("H").minus(tru.get("H")).print(5,5);
		
		System.out.println("P's");
		emp.get("pinv").print(5, 5);
		tru.get("pinv").print(5, 5);
		emp.get("pinv").minus(tru.get("pinv")).print(5, 5);
		
		System.out.println("SVDs");
		emp.get("s_values").print(5, 5);
		tru.get("s_values").print(5, 5);
		
		
		System.out.println("U's");
		emp.get("U").print(5, 5);
		tru.get("U").print(5, 5);
		
		System.out.println("VT's");

		emp.get("VT").print(5, 5);
		tru.get("VT").print(5, 5);
		
		System.out.println("Hsigma=1 error");
		emp.get("1").minus(tru.get("1")).print(5,5);
	}

}
