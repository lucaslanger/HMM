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
	private int maxQuery;
	private int maxExp;
	private int maxStates;

	public static void main(String[] args){
		int hSize = 80;
		int basisSize = 40;
 	
		Analysis a = new Analysis(hSize, basisSize);
	}
	
	public Analysis(int hSize, int basisSize){
		this.hSize = hSize;
		this.basisSize = basisSize;
				
		int firstLoop = 18;
		int secondLoop = 12;
		this.maxStates = firstLoop + secondLoop - 1;
		
		double selfTransition = 0.05;
		this.h = this.makeLabyrinth(firstLoop,secondLoop,selfTransition);
		this.tru = h.singledataSpectralTrue(hSize, basisSize, this.maxStates);
		this.maxExp = (int) this.tru.get("max").get(0, 0);
		this.maxQuery = (int) Math.pow(2,this.maxExp);
		
		/*
		int rep1 = 1;
		int amountOfData = 70;
		this.fixedSizePlots(rep1, amountOfData);
		
		int rep2 = 50;
		this.plotBaseDifferences( h, hSize, 40 , rep2);
		*/
		
		int rep3 = 15;
		this.sizeOfModelPlots(rep3);
		
		/*
		int[] dataAmount = new int[]{50,70,100,150,200,500,1000,10000};
		for (int i = 0; i < dataAmount.length; i++) {
			for (int j = 0; j < repeats; j++) {
				//Take out redundencies if possible
			}
		}
		*/
		
	}
	
	public void fixedSizePlots(int trials, int amountOfData){		
		for (int i = 0; i < trials; i++) {
			this.empArray.add(h.singledataSpectralEmperical(hSize, amountOfData, basisSize, this.maxStates));
		}
		
		this.conditionalPlots(trials,1,50);
		//System.out.println("Conditional Done");
		this.compareSigmaError();
		//System.out.println("Power Sigma Done");
		this.compareH_Hbar(h, 1);
		//System.out.println("H Hbar Done");
		this.compareASigmas();
		//System.out.println("SigmaError Done");
		this.compareQueryErrors();
		//System.out.println("Query Done");
	}
	
	public void plotBaseDifferences( HMM hmm, int hankelSize, int basisSize, int repeats){
		int[] dataAmount = new int[]{50,70,100,150,200,500,1000,10000};
		
		HashMap<String, Matrix> d = hmm.singledataSpectralTrue(hankelSize, basisSize, this.maxStates);
		HashMap<String, Matrix> emp;
		
		int maxExp = (int) d.get("max").get(0,0);
		double[][] dataSize = new double[maxExp+1][dataAmount.length];
		double[][] errors = new double[maxExp+1][dataAmount.length];
		
		double error = 0, empProb, truProb;
		int maxQuery = hankelSize;
		int exp;
		
		for (int j = 0; j < dataAmount.length; j++){
			for (int z = 0; z < repeats; z++){
				emp = hmm.singledataSpectralEmperical(hankelSize, dataAmount[j], basisSize, this.maxStates);	
				for (int i = 0; i <= maxExp; i++){
					exp = (int) Math.pow(2, i);
					dataSize[i][j] = dataAmount[j];
					for (int j2 = 0; j2 < maxQuery; j2++){
						truProb = HelperFunctions.probabilityQuery(d, d.get("a0"), d.get("ainf"), j2, exp, 2, true);
						empProb = HelperFunctions.probabilityQuery(emp, emp.get("a0"), emp.get("ainf"), j2, exp, 2, true);
						error = computeError(truProb, empProb);
						errors[i][j] += error;
					}
				}
			}
		}
		
		for (int i = 0; i <= maxExp; i++) {
			for (int j = 0; j < dataAmount.length; j++) {
				errors[i][j] /= repeats;
			}
		}
		
		Matrix visualErrors = new Matrix(errors);
		visualErrors.print(5, 5);
		
		HelperFunctions.outputData(pltFolder + "BaseComp_Area", "X:#Data Seen Y:Fnorm","", dataSize, errors );
	}
	
	private double computeError(double truProb, double empProb) {
		double error =  Math.abs( truProb - empProb );
		
		return error;
	}

	public void conditionalPlots(int trials, int traj, int maxAhead){
		
		double[][] xaxis = new double[traj][maxAhead];
		double[][] queryArrayTru = new double[traj][maxAhead];
		for (int i = 0; i < traj; i=i+1) {
			queryArrayTru[i] = conditionalQuery(this.tru,i,maxAhead);
			xaxis[i] = HelperFunctions.incArray(maxAhead);
			//System.out.println(HelperFunctions.sumArray(queryArray[i]));  
		}
		Matrix truPredictions = new Matrix(queryArrayTru);

		double[][] queryArrayEmp = new double[traj][maxAhead];
		Matrix queryEmpAvg = null, qE = null;
		double[][] errorArray;
		Matrix error, errorAbs, avgError = null;
		
		for (int i = 0; i < trials; i++) {
			for (int j = 0; j < traj; j++) {
				queryArrayEmp[j] = conditionalQuery(this.empArray.get(i), j, maxAhead);
			}
			qE = new Matrix(queryArrayEmp);
			if (queryEmpAvg != null){
				queryEmpAvg = queryEmpAvg.plus(qE);
			}
			else{
				queryEmpAvg = qE;
			}
			
			error = qE.minus(truPredictions);
			
			errorArray = new double[traj][maxAhead];
			for (int j = 0; j < traj; j++) {
				for (int j2 = 0; j2 < maxAhead; j2++) {
					errorArray[j][j2] = Math.abs(error.get(j,j2));
				}
			}
			errorAbs = new Matrix(errorArray);
			
			if (avgError != null){
				avgError = avgError.plus(errorAbs);
			}
			else{
				avgError = errorAbs;
			}
			
		}
		
		avgError = avgError.times(1.0/trials);
		queryEmpAvg = queryEmpAvg.times(1.0/trials);
		
		HelperFunctions.outputData(pltFolder + "ConditionalError", "x:Traj Length y:|f_k(x)-fhat_k(x)|", "", xaxis, avgError.getArrayCopy());
		HelperFunctions.outputData(pltFolder + "ConditionalEmp", "x:Traj Length y:fhat_k(x)", "", xaxis, queryEmpAvg.getArrayCopy());
		HelperFunctions.outputData(pltFolder + "ConditionalTrue", "x:Traj Length y:f_k(x)", "", xaxis, truPredictions.getArrayCopy());
	}
	
	public void compareH_Hbar(HMM hmm, int repeats){
		
		int[] sizes = new int[]{100,500,1000,2000,5000,10000};
		
		HashMap<String, Matrix> emp;
		
		Matrix H;
		Matrix Hbar;
		
		int trialsize = sizes.length;
		
		double[][] dataSize = new double[1][trialsize];
		double[][] error = new double[1][trialsize];
		double avgError, e;
	
		for (int i = 0; i < trialsize; i++) {
			avgError = 0;
			for (int j = 0; j < repeats; j++) {
				emp = hmm.singledataSpectralEmperical(hSize, sizes[i], basisSize, this.maxStates);
				
				H = this.tru.get("H");
				Hbar = emp.get("H");
				
				e = H.minus(Hbar).normF();
				avgError += e;
			}
			avgError /= repeats;
			dataSize[0][i] = sizes[i];
			error[0][i] = avgError;
		}
			
		HelperFunctions.outputData(pltFolder + "True_H_vs_Emp", "X:#Data Seen Y:Fnorm","", dataSize, error );
		
	}
	
	public void compareASigmas(){
		HashMap<String, Matrix> emp;
		
		int maxQuery = (int) (this.tru.get("max").get(0,0) ); 
		
		double[][] sigmaNumber = new double[1][maxQuery];
		double[][] errors = new double[1][maxQuery];
		
		Matrix h_sigma_true, h_sigma_exp, r;
		int pow;
		
		for (int j = 0; j < this.empArray.size(); j++) {
			emp = this.empArray.get(j);
			for (int i = 0; i < maxQuery; i++) {
				pow = (int) Math.pow(2, i);
				h_sigma_true = this.tru.get( Integer.toString(pow) );
				h_sigma_exp = emp.get( Integer.toString(pow) );
				r = h_sigma_true.minus( h_sigma_exp );
				errors[0][i] += r.normF();
			}
		}
		
		for (int i = 0; i < maxQuery; i++) {
			pow = (int) Math.pow(2, i);
			h_sigma_true = this.tru.get( Integer.toString(pow) );
			sigmaNumber[0][i] = pow;
			errors[0][i] /= (this.empArray.size()*h_sigma_true.normF());
		}
		
		HelperFunctions.outputData(pltFolder + "True_Ax_vs_Emp", "X:Sigma Y:(T_Ax-E_Ax).Fnorm/T_Ax.Fnorm","", sigmaNumber, errors );
		// Add file containing error analysis for alphaInf and alpha0?
	}
	
	public void compareSigmaError(){
		int m = (int) (this.tru.get("max").get(0,0) - 1); 
		
		double[][] sigmaNumber = new double[1][m];
		double[][] errors = new double[1][m];
		
		Matrix temp1, temp2, r;
		int pow;
		
		for (int j = 0; j < this.empArray.size(); j++) {
				
			for (int i = 0; i < m; i++) {
				pow = (int) Math.pow(2, i);
				temp1 = this.empArray.get(j).get( Integer.toString( pow ) );
				temp1 = HelperFunctions.matrixPower( temp1 , 2);
				temp2 = this.empArray.get(j).get( Integer.toString(pow*2) );
				r = temp2.minus( temp1 ) ;	
					
				errors[0][i] += r.normF();
			}
			
		}
		
		Matrix h_sigma_true;
		for (int i = 0; i < m; i++) {
			pow = (int) Math.pow(2, i+1);
			h_sigma_true = this.tru.get( Integer.toString(pow) );
			sigmaNumber[0][i] = pow;
			errors[0][i] /= (this.empArray.size()*h_sigma_true.normF()*hSize*hSize);
		}
		
		HelperFunctions.outputData(pltFolder + "(Ax)^2_v.s A(x^2)", "X:Sigma Y:(T_Ax-E_Ax).Fnorm/T_Ax.Fnorm","", sigmaNumber, errors );

	}
	
	public void compareQueryErrors(){
		int maxexp = (int) this.tru.get("max").get(0,0); //same for tru by construction
		int maxpow = (int) Math.pow(2, maxexp);
		
		int extraQuery = 50;
		int maxquery = maxpow + extraQuery;
		
		double[][] queries = new double[9][maxquery];
		double[][] errors = new double[9][maxquery];
		
		double[][] baseQueries = new double[maxexp][maxquery];
		double[][] x_base_Queries = new double[maxexp][maxquery];
		
		Matrix a0emp, ainfemp, empQF, empQB;
		Matrix a0tru, ainftru;
		double truProbQF, truProbQB , empProbQF, empProbQB, empProbP;
		
		HashMap<String, Matrix> emp;
		for (int i = 0; i < maxquery ; i++) {
			//truF = HelperFunctions.matrixQuery(tru, i, this.maxPower, 2, true);
			//truB = HelperFunctions.matrixQuery(tru, i, this.maxPower, 2, false);
			a0tru = tru.get("a0");
			ainftru = tru.get("ainf");
	
			truProbQF = HelperFunctions.probabilityQuery(tru, a0tru, ainftru, i, maxpow, 2, true);
			truProbQB = HelperFunctions.probabilityQuery(tru, a0tru, ainftru, i, maxpow, 2, false);
			//truProbF.minus(truProbB).print(5, 5); //Always 0 which makes sense
			
			for (int j = 0; j < this.empArray.size(); j++) {	
				emp = this.empArray.get(j);
		
				empQF = HelperFunctions.matrixQuery(emp, i, maxpow, 2, true);
				empQB = HelperFunctions.matrixQuery(emp, i, maxpow, 2, false);
				//empP = HelperFunctions.matrixPower(emp.get("1"), i);										//inefficient, if slow optimize later
				
				a0emp = emp.get("a0");
				ainfemp = emp.get("ainf");
				
				/*empProbQF = a0emp.times(empQF).times(ainfemp);
				empProbQB = a0emp.times(empQB).times(ainfemp);
				empProbP = a0emp.times(empP).times(ainfemp);
				*/
				
				empProbQF = HelperFunctions.probabilityQuery(emp, a0emp, ainfemp, i, maxpow, 2, true);
				empProbQB = HelperFunctions.probabilityQuery(emp, a0emp, ainfemp, i, maxpow, 2, false);
				empProbP = HelperFunctions.probabilityQuery(emp, a0emp, ainfemp, i, 1, 2, true);
				
				errors[0][i] += Math.abs(truProbQF - empProbQF);	//Tru v.s Base 
				errors[1][i] += truProbQF - empProbQF;

				errors[2][i] += Math.abs(truProbQF - empProbP);		//Tru v.s Naive 
				errors[3][i] += truProbQF - empProbP ;
								
				errors[4][i] += Math.abs(empProbQF - empProbQB );	// Comm Error
				errors[5][i] += empProbQF - empProbQB ;
				
				errors[6][i] += empQF.minus(empQB).normF();			// Matrix Comm Error
				
				/*
				double r7 = Math.max( Math.abs(truProbQF), Math.abs(empProbQF) );
				double r8 = Math.max( Math.abs(truProbQF), Math.abs(empProbP) );
				
				errors[7][i] += Math.abs(truProbQF - empProbQF);	
				errors[8][i] += Math.abs(truProbQF - empProbP);
				*/
				
				double pq;
				for (int k = 0; k < baseQueries.length; k++) {
					pq = Math.abs(HelperFunctions.probabilityQuery(emp, a0emp, ainfemp, i, (int) Math.pow(2,k), 2, true) - truProbQF);
					baseQueries[k][i] += pq;
				}
			}
			
			for (int j = 0; j < x_base_Queries.length; j++){
				baseQueries[j][i] /= (this.empArray.size() );
				x_base_Queries[j] = HelperFunctions.incArray(maxquery);
			}
			
			for (int c = 0; c < 9; c++) {
				errors[c][i] /= (this.empArray.size());
				queries[c][i] = i;
			}
		}
		
		HelperFunctions.outputData(pltFolder + "Query_Errors_Base", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,0,2), Arrays.copyOfRange(errors,0,2) );
		HelperFunctions.outputData(pltFolder + "Query_Errors_Naive", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,2,4), Arrays.copyOfRange(errors,2,4) );
		HelperFunctions.outputData(pltFolder + "Comm_Query_Error", "X:Sigma Y:a0(A16A1-A1A16)aI","", Arrays.copyOfRange(queries,4,6), Arrays.copyOfRange(errors,4,6) );
		HelperFunctions.outputData(pltFolder + "Comm_Matrix_Error", "X:Sigma Y:(A16A1-A1A16).Fnorm","", Arrays.copyOfRange(queries,6,7), Arrays.copyOfRange(errors,6,7) );
		HelperFunctions.outputData(pltFolder + "Base_Errors","X:Sigma Y: Error" ,"", x_base_Queries, baseQueries);
		
		double[][] ebase = Arrays.copyOfRange(errors,0, 1);
		double[][] enaive = Arrays.copyOfRange(errors, 2, 3);
		double[][] ejoint = new double[][]{ebase[0], enaive[0]};
		
		double[][] qbase = Arrays.copyOfRange(queries, 2, 3);
		double[][] qnaive = Arrays.copyOfRange(queries, 2, 3);
		double[][] qjoint = new double[][]{qbase[0], qnaive[0]};
		HelperFunctions.outputData(pltFolder + "QError_Base_vs_Naive", "X:Sigma Y:|f(x)-fhat(x)|","",qjoint,ejoint  );
		
		System.out.println("Highest Base = " + Integer.toString(maxpow));
		System.out.println( HelperFunctions.sumArray(errors[0]) );
		System.out.println("Naive Max-Base = 1");
		System.out.println( HelperFunctions.sumArray(errors[2]) );

		
	}
	
	public HMM makeHMM(){
		double[][] p = { {1}, {0}};
		double[][] t = { {0.5,0.45}, {0.3,0.67} };
		double[][] o = { {0,1}, {0,1} };
		double[][] e = { {0.05}, {0.03} };
		
		Matrix T = new Matrix( t );
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		Matrix E = new Matrix( e );
		
		HMM h = new HMM(T,O,P,E);	
		return h;
	}
	
	public HMM makeLabyrinth(int loop1, int loop2 , double selfTransitionP){
		int states = loop1 + loop2 - 1;
		HashMap<Integer, Double> termStates = new HashMap<Integer, Double>();
		int door1 = 0;
		int door2 = loop1/2 + loop2/2;
		termStates.put(door1, .6);
		termStates.put(door2, .4);
		
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
		
		Matrix T = new Matrix( t ).transpose();
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p );
		Matrix E = new Matrix( e );
		
		HMM l = new HMM(T, O, P, E);
		
		return l;
		
	}
	
	public double[] conditionalQuery(HashMap<String, Matrix> learned, int k, int maxAhead){
		int maxpow = (int) Math.pow(2,learned.get("max").get(0, 0));
		Matrix alpha_0 = learned.get("a0");
		Matrix alpha_inf = learned.get("ainf");
		//Matrix Ak = HelperFunctions.matrixQuery(learned, k, 2, true);
		Matrix alpha_k = HelperFunctions.alphaKQuery(learned, alpha_0, k, maxpow, 2);//alpha_0.times(Ak);
		
		int nstates = alpha_0.getArray()[0].length;
		
		Matrix mid = Matrix.identity(nstates, nstates).minus(learned.get("1") );
		double normalizer = alpha_k.times( mid.inverse() ).times( alpha_inf ).get(0,0);
		double jointProb;
		
		double[] pA = new double[maxAhead];
		int maxbase = 1;
		for (int i = 0; i < pA.length; i++) {
			//Aquery = HelperFunctions.matrixQuery(learned, i, 2,true);				//Inefficient way to multiply fix to optimize if slow
			//jointProb = alpha_k.times(Aquery).times(alpha_inf).get(0,0);
			jointProb = HelperFunctions.probabilityQuery(learned, alpha_k, alpha_inf, i, maxbase ,2, true);
			pA[i] = jointProb/normalizer;			
		}
		
		return pA;
	}
	
	public void sizeOfModelPlots(int repeats){
		//HashMap<Integer, ArrayList<ArrayList<HashMap<String, Matrix> >>> data = null;
		
		int[] dataAmount = new int[]{50,70,100,150,200,500,1000,2000,4000};
		double[][] plotErrors = new double[this.maxExp+1][dataAmount.length];
		double[][] plotArgForErrors = new double[this.maxExp+1][dataAmount.length];
		
		double[][] xaxis = new double[this.maxExp+1][dataAmount.length];
		for (int i = 0; i < xaxis.length; i++) {
			for (int j = 0; j < dataAmount.length; j++) {
				xaxis[i][j] = dataAmount[j];
			}
		}
		
		int baseSize;
		for (int c = 0; c <= maxExp; c++) {

			baseSize = (int) Math.pow(2, c);;
			System.out.print("Base ");
			System.out.println(baseSize);
			
			double[][] errors = new double[dataAmount.length][this.maxStates];
			double truQuery, empQuery, error;
			
			HashMap<String, Matrix> t;
			for (int i = 0; i < dataAmount.length; i++){
				for (int j = 1; j <= errors[0].length; j++){	
					for (int z = 0; z < repeats; z++){
						t = this.h.singledataSpectralEmperical(this.hSize, dataAmount[i], this.basisSize, j);
						for (int q = 0; q < this.maxQuery; q++){
							empQuery = HelperFunctions.probabilityQuery(t, t.get("a0"),  t.get("ainf"), q, baseSize, 2, true);
							truQuery = HelperFunctions.probabilityQuery(this.tru, this.tru.get("a0"),  this.tru.get("ainf"), q, baseSize, 2, true);
							error = computeError(truQuery, empQuery);
							errors[i][j-1] += error;
						}
					}
					errors[i][j-1] /= repeats;
				}
			}
			
			//Matrix printError = new Matrix(errors);
			//printError.print(5, 5);
			
			double[] argMinArray = new double[dataAmount.length];
			double[] errorMinArray = new double[dataAmount.length];
			
			for (int i = 0; i < errorMinArray.length; i++) {
				argMinArray[i] =  HelperFunctions.getMinValue( errors[i] );
				errorMinArray[i] = HelperFunctions.getArgMin( errors[i] ); 
				/*System.out.println("Min Errors");
				System.out.println(errorMinArray[i]);
				System.out.println(argMinArray[i]);
				*/
			}
			
			plotErrors[c] = errorMinArray;
			plotArgForErrors[c] = argMinArray;
			
		}
		
		Matrix printBestBaseErrors = new Matrix(plotErrors);
		Matrix printBestBaseArg = new Matrix(plotArgForErrors);
		
		printBestBaseErrors.print(5, 5);
		printBestBaseArg.print(5, 5);
	
		HelperFunctions.outputData(pltFolder + "MinError_Dif_Bases", "X: Data, Y:Min_over_#states", "", xaxis, plotErrors);
		HelperFunctions.outputData(pltFolder + "ArgMin_Dif_Bases", "X: Data, Y:ArgMin_over_#states", "", xaxis, plotArgForErrors);
	}

	
	public void debugHComparisons(HashMap<String, Matrix> machine ){
		System.out.println("H error");
		machine.get("H").minus(tru.get("H")).print(5,5);
		
		System.out.println("P's");
		machine.get("pinv").print(5, 5);
		tru.get("pinv").print(5, 5);
		machine.get("pinv").minus(tru.get("pinv")).print(5, 5);
		
		System.out.println("SVDs");
		machine.get("s_values").print(5, 5);
		tru.get("s_values").print(5, 5);	
		
		System.out.println("U's");
		machine.get("U").print(5, 5);
		tru.get("U").print(5, 5);
		
		System.out.println("VT's");

		machine.get("VT").print(5, 5);
		tru.get("VT").print(5, 5);
		
		System.out.println("Asigma=1 error");
		machine.get("1").minus(tru.get("1")).print(5,5);
	}

}
