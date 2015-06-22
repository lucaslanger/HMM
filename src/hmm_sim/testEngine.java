package hmm_sim;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import Jama.Matrix;

public class testEngine{
	
	public String pltFolder;
	private String fileNameOfDataSet;

	private HashMap<Integer, HankelSVDModel[]> empericalModels;
	private QueryEngine[][] fixedSizeQueryEngines;
	private QueryEngine[][][] anySizeQueryEngines;
	
	private HankelSVDModel trueModel;
	private QueryEngine trueQueryEngine;

	private int base;
	private int basisSize;
	private int maxQuery;
	private int maxStates;

	
	
	public static void main(String[] args){
 	
		testEngine a = new testEngine();
	}
	
	public testEngine(String fileNameOfDataSet, int basisSize, int base, int numberPerTrajectorySize){
		this.fileNameOfDataSet = fileNameOfDataSet;
		this.basisSize = basisSize;
		this.base = base;
		
		this.trueModel = this.readTrueModel(fileNameOfDataSet);
		this.maxStates = this.trueModel.getRank();
		this.trueQueryEngine = this.trueModel.buildHankelBasedModel(basisSize, base, this.maxStates);
		
		this.empericalModels = this.readEmpericalModels(fileNameOfDataSet, numberPerTrajectorySize);
				
		this.anySizeQueryEngines = this.getAllSizesQueryEngines(100);
		
		this.fixedSizePlots();
		System.out.println("Done Fixed Plots");
		
		this.plotBaseDifferences( );
		System.out.println("Done Base Differences");
		
		this.sizeOfModelPlots( );
		System.out.println("Done Model Differences");
		

	}
	
	private HashMap<Integer, QueryEngine[][]> getAllSizesQueryEngines(int numberOfTrajectoriesFromEachSize){
		String[] files;
		HankelSVDModel h;
		HashMap<Integer, QueryEngine[][]> dataSizeToModels = null;
		QueryEngine[][] Q;
		QueryEngine q;

		try{
			for (String file: files) {
				int trajectoryLength = Integer.parseInt(file);
				ObjectInputStream ois = new ObjectInputStream( new FileInputStream(file) );
				for (int i = 0; i < numberOfTrajectoriesFromEachSize; i++) {
					h = (HankelSVDModel) ois.readObject();
					for (int j = 0; j < this.maxStates; j++) {
						q = h.buildHankelBasedModel(i, numberOfTrajectoriesFromEachSize, j);
						Q[i][j] = q;
					}
				}
				dataSizeToModels.put(trajectoryLength, Q);
			}
			return dataSizeToModels;
		}
		catch(Exception e){
			e.printStackTrace();
			return null;
		}
		
	}

	public HankelSVDModel readTrueModel(String f){
		HankelSVDModel t;
		//Search for file with True model
		ObjectInputStream ois;
		try {
			ois = new ObjectInputStream(new FileInputStream(f));
			t = (HankelSVDModel) ois.readObject();
			return t;
		} catch (IOException e) {
			System.out.println("Problem fetching true Model");
			e.printStackTrace();
			return null;
		} catch (ClassNotFoundException e) {
			System.out.println("Couldn't find the file");
			e.printStackTrace();
			return null;
		}
		
	}
	
	public HashMap<Integer, HankelSVDModel[]> readEmpericalModels(String f, int numberPerTrajectorySize){
		String[] files;
		HashMap<Integer, HankelSVDModel[]> h;
		
		HankelSVDModel[] empericalModelsFixedTrajectories;
		for (String file: files) {	//
			int numberOfTrajectories = Integer.toString(file);
			ObjectInputStream ois = new ObjectInputStream(new FileInputStream(file));
			
			empericalModelsFixedTrajectories = new HankelSVDModel[numberPerTrajectorySize];
			for (int i = 0; i < empericalModelsFixedTrajectories.length; i++) {
				empericalModelsFixedTrajectories[i] = (HankelSVDModel) ois.readObject();
			}
			h.put(numberOfTrajectories, empericalModelsFixedTrajectories);
		}
		return h;
		
	}
	
	public void fixedSizePlots(){		
		
		int maxAhead = 50;
		this.conditionalPlots(maxAhead);
		
		this.compareSigmaError();
		
		this.compareH_Hbar(5);
		
		if(this.maxStates == nStates){
			this.compareASigmas();
		}
		
		System.out.println("SigmaError Done");
		this.compareQueryErrors();
		System.out.println("Query Done");
		
		int num_lines = 10;
		int amountOfData2 = 100;
		this.plotTrialsModelSize(num_lines, amountOfData2);
		System.out.println("Done Multiple Lines Model Error");
	}
	
	public void plotBaseDifferences( int hankelSize, int basisSize, int repeats){
		
		HashMap<String, Matrix> emp;
		
		double[][] dataSize = new double[this.maxExp+1][this.dataSizes.length];
		double[][] errors = new double[this.maxExp+1][this.dataSizes.length];
		
		double error = 0, empProb, truProb;
		int maxQuery = hankelSize;
		int exp;
		
		for (int j = 0; j < this.dataSizes.length; j++){
			for (int z = 0; z < repeats; z++){
				emp = this.h.singledataSpectralEmperical(hankelSize, this.dataSizes[j], basisSize, this.maxStates);
				for (int i = 0; i <= this.maxExp; i++){
					exp = (int) Math.pow(2, i);
					dataSize[i][j] = this.dataSizes[j];
					for (int j2 = 0; j2 < maxQuery; j2++){
						truProb = Analysis.probabilityQuery(this.tru, this.tru.get("a0"), this.tru.get("ainf"), j2, exp, 2, true);
						empProb = Analysis.probabilityQuery(emp, emp.get("a0"), emp.get("ainf"), j2, exp, 2, true);
						error = computeError(truProb, empProb);
						errors[i][j] += error;
					}
				}
			}
		}
		
		for (int i = 0; i <= this.maxExp; i++) {
			for (int j = 0; j < this.dataSizes.length; j++) {
				errors[i][j] /= repeats;
			}
		}
		Analysis.outputData(pltFolder + "BaseComp_Area", "X:#Data Seen Y:Fnorm","", dataSize, errors );

		
		System.out.println("");
		System.out.println("Base Comp Errors Modelsize=" + Integer.toString(this.maxStates));
		System.out.println("Downwards: BASE, SideWays: #DATA");
		Matrix visualErrors = new Matrix(errors);
		visualErrors.print(5, 5);
		
	}
	
	private double computeError(double truProb, double empProb) {
		double error =  Math.abs( truProb - empProb );
		return error;
	}

	public void conditionalPlots(int maxAhead){
		
		double[][] xaxis = new double[traj][maxAhead];
		double[][] queryArrayTru = new double[traj][maxAhead];
		for (int i = 0; i < traj; i=i+1) {
			queryArrayTru[i] = conditionalQuery(this.tru,i,maxAhead);
			xaxis[i] = Analysis.incArray(maxAhead);
			//System.out.println(Analysis.sumArray(queryArray[i]));  
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
		
		testEngine.outputData(pltFolder + "ConditionalError", "x:Traj Length y:|f_k(x)-fhat_k(x)|", "", xaxis, avgError.getArrayCopy());
		testEngine.outputData(pltFolder + "ConditionalEmp", "x:Traj Length y:fhat_k(x)", "", xaxis, queryEmpAvg.getArrayCopy());
		testEngine.outputData(pltFolder + "ConditionalTrue", "x:Traj Length y:f_k(x)", "", xaxis, truPredictions.getArrayCopy());
	}
	
	public void compareH_Hbar(int repeats){
		
		HashMap<String, Matrix> emp;
		
		Matrix H = this.tru.get("H");
		Matrix Hbar;
		
		double[][] dataSize = new double[1][this.dataSizes.length];
		double[][] error = new double[1][this.dataSizes.length];
		double avgError, e;
	
		for (int i = 0; i < this.dataSizes.length; i++) {
			avgError = 0;
			for (int j = 0; j < repeats; j++) {
				emp = this.h.singledataSpectralEmperical(hSize, this.dataSizes[i], basisSize, this.maxStates);
				Hbar = emp.get("H");
				
				e = H.minus(Hbar).normF();
				avgError += e;
			}
			avgError /= repeats;
			dataSize[0][i] = this.dataSizes[i];
			error[0][i] = avgError;
		}
			
		Analysis.outputData(pltFolder + "True_H_vs_Emp", "X:#Data Seen Y:Fnorm","", dataSize, error );
		
	}
	
	public void compareASigmas(){
		HashMap<String, Matrix> emp;
		
		double[][] sigmaNumber = new double[1][this.maxExp];
		double[][] errors = new double[1][this.maxExp];
		
		Matrix a_sigma_true, a_sigma_exp, r;
		int pow;
		
		for (int j = 0; j < this.empArray.size(); j++) {
			emp = this.empArray.get(j);
			for (int i = 0; i < errors[0].length; i++) {
				pow = (int) Math.pow(2, i);
				a_sigma_true = this.tru.get( Integer.toString(pow) );
				a_sigma_exp = emp.get( Integer.toString(pow) );

				r = a_sigma_true.minus( a_sigma_exp );
				errors[0][i] += r.normF();
			}
		}
		
		for (int i = 0; i < errors[0].length; i++) {
			pow = (int) Math.pow(2, i);
			a_sigma_true = this.tru.get( Integer.toString(pow) );
			sigmaNumber[0][i] = pow;
			errors[0][i] /= (this.empArray.size()*a_sigma_true.normF());
		}
		
		Analysis.outputData(pltFolder + "True_Ax_vs_Emp", "X:Sigma Y:(T_Ax-E_Ax).Fnorm/T_Ax.Fnorm","", sigmaNumber, errors );
		// Add file containing error analysis for alphaInf and alpha0?
	}
	
	public void compareSigmaError(){
		
		int maxExpSquareComparison = this.maxExp-1;
		
		double[][] sigmaNumber = new double[1][maxExpSquareComparison];
		double[][] errors = new double[1][maxExpSquareComparison];
		
		Matrix temp1, temp2, r;
		int pow;
		
		for (int j = 0; j < this.empArray.size(); j++) {
				
			for (int i = 0; i < maxExpSquareComparison; i++) {
				pow = (int) Math.pow(2, i);
				temp1 = this.empArray.get(j).get( Integer.toString( pow ) );
				temp1 = Analysis.matrixPower( temp1 , 2);
				temp2 = this.empArray.get(j).get( Integer.toString(pow*2) );
				r = temp2.minus( temp1 ) ;	
					
				errors[0][i] += r.normF();
			}
			
		}
		
		Matrix h_sigma_true;
		for (int i = 0; i < maxExpSquareComparison; i++) {
			pow = (int) Math.pow(2, i+1);
			h_sigma_true = this.tru.get( Integer.toString(pow) );
			sigmaNumber[0][i] = pow;
			errors[0][i] /= (this.empArray.size()*h_sigma_true.normF()*hSize*hSize);
		}
		
		Analysis.outputData(pltFolder + "(Ax)^2_v.s A(x^2)", "X:Sigma Y:(T_Ax-E_Ax).Fnorm/T_Ax.Fnorm","", sigmaNumber, errors );

	}
	
	public void compareQueryErrors(){		
		
		double[][] queries = new double[9][this.maxQuery];
		double[][] errors = new double[9][this.maxQuery];
		
		double[][] baseQueries = new double[this.maxExp][this.maxQuery];
		double[][] x_base_Queries = new double[this.maxExp][this.maxQuery];
		
		Matrix a0emp, ainfemp, empQF, empQB;
		Matrix a0tru, ainftru;
		double truProbQF, truProbQB, truProbP , empProbQF, empProbQB, empProbP;
		
		HashMap<String, Matrix> emp;
		for (int i = 0; i < this.maxQuery ; i++) {
	
			a0tru = tru.get("a0");
			ainftru = tru.get("ainf");
	
			truProbQF = Analysis.probabilityQuery(tru, a0tru, ainftru, i, this.maxPower, 2, true);
			truProbQB = Analysis.probabilityQuery(tru, a0tru, ainftru, i, this.maxPower, 2, false);
			truProbP = Analysis.probabilityQuery(tru, a0tru, ainftru, i, 1, 2, false);
			//System.out.println(truProbQF - truProbQB);//Always 0 which makes sense
			//System.out.println(truProbP - truProbQF);
			
			for (int j = 0; j < this.empArray.size(); j++) {	
				emp = this.empArray.get(j);
		
				empQF = Analysis.matrixQuery(emp, i, this.maxPower, 2, true);
				empQB = Analysis.matrixQuery(emp, i, this.maxPower, 2, false);
				//empP = Analysis.matrixPower(emp.get("1"), i);										//inefficient, if slow optimize later
				
				a0emp = emp.get("a0");
				ainfemp = emp.get("ainf");
				
				/*empProbQF = a0emp.times(empQF).times(ainfemp);
				empProbQB = a0emp.times(empQB).times(ainfemp);
				empProbP = a0emp.times(empP).times(ainfemp);
				*/
				
				empProbQF = Analysis.probabilityQuery(emp, a0emp, ainfemp, i, this.maxPower, 2, true);
				empProbQB = Analysis.probabilityQuery(emp, a0emp, ainfemp, i, this.maxPower, 2, false);
				empProbP = Analysis.probabilityQuery(emp, a0emp, ainfemp, i, 1, 2, true);
				
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
					pq = Math.abs(Analysis.probabilityQuery(emp, a0emp, ainfemp, i, (int) Math.pow(2,k), 2, true) - truProbQF);
					baseQueries[k][i] += pq;
				}
			}
			
			for (int j = 0; j < x_base_Queries.length; j++){
				baseQueries[j][i] /= (this.empArray.size() );
				x_base_Queries[j] = Analysis.incArray(this.maxQuery);
			}
			
			for (int c = 0; c < 9; c++) {
				errors[c][i] /= (this.empArray.size());
				queries[c][i] = i;
			}
		}
		
		Analysis.outputData(pltFolder + "Query_Errors_Base", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,0,2), Arrays.copyOfRange(errors,0,2) );
		Analysis.outputData(pltFolder + "Query_Errors_Naive", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,2,4), Arrays.copyOfRange(errors,2,4) );
		Analysis.outputData(pltFolder + "Comm_Query_Error", "X:Sigma Y:a0(A16A1-A1A16)aI","", Arrays.copyOfRange(queries,4,6), Arrays.copyOfRange(errors,4,6) );
		Analysis.outputData(pltFolder + "Comm_Matrix_Error", "X:Sigma Y:(A16A1-A1A16).Fnorm","", Arrays.copyOfRange(queries,6,7), Arrays.copyOfRange(errors,6,7) );
		Analysis.outputData(pltFolder + "Base_Errors","X:Sigma Y: Error" ,"", x_base_Queries, baseQueries);
		
		double[][] ebase = Arrays.copyOfRange(errors,0, 1);
		double[][] enaive = Arrays.copyOfRange(errors, 2, 3);
		double[][] ejoint = new double[][]{ebase[0], enaive[0]};
		
		double[][] qbase = Arrays.copyOfRange(queries, 2, 3);
		double[][] qnaive = Arrays.copyOfRange(queries, 2, 3);
		double[][] qjoint = new double[][]{qbase[0], qnaive[0]};
		Analysis.outputData(pltFolder + "QError_Base_vs_Naive", "X:Sigma Y:|f(x)-fhat(x)|","",qjoint,ejoint  );
		
		/*	//Print out area under curve between naive and base method
		System.out.println("Highest Max-Base = " + Integer.toString(this.maxPower));
		System.out.println( Analysis.sumArray(errors[0]) );
		System.out.println("Naive Max-Base = 1");
		System.out.println( Analysis.sumArray(errors[2]) );
		*/
		
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
		boolean debug = false;
		
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
		
		Matrix T = new Matrix( t ).transpose();
		
		Matrix O = new Matrix( o );
		Matrix P = new Matrix( p ).transpose();
		Matrix E = new Matrix( e );
		
		HMM l = new HMM(T, O, P, E);
		
		return l;
		
	}
	
	public double[] conditionalQuery(HashMap<String, Matrix> learned, int k, int maxAhead){
		int maxpow = (int) Math.pow(2,learned.get("max").get(0, 0));
		Matrix alpha_0 = learned.get("a0");
		Matrix alpha_inf = learned.get("ainf");
		//Matrix Ak = Analysis.matrixQuery(learned, k, 2, true);
		Matrix alpha_k = Analysis.alphaKQuery(learned, alpha_0, k, maxpow, 2);//alpha_0.times(Ak);
		
		int nstates = alpha_0.getArray()[0].length;
		
		Matrix mid = Matrix.identity(nstates, nstates).minus(learned.get("1") );
		double normalizer = alpha_k.times( mid.inverse() ).times( alpha_inf ).get(0,0);
		double jointProb;
		
		double[] pA = new double[maxAhead];
		int maxbase = 1;
		for (int i = 0; i < pA.length; i++) {
			jointProb = Analysis.probabilityQuery(learned, alpha_k, alpha_inf, i, maxbase ,2, true);
			pA[i] = jointProb/normalizer;			
		}
		
		return pA;
	}
	
	public void sizeOfModelPlots(int repeats, boolean debug){
		//HashMap<Integer, ArrayList<ArrayList<HashMap<String, Matrix> >>> data = null;
		
		double[][] plotErrors = new double[this.maxExp+1][this.dataSizes.length];
		double[][] plotArgForErrors = new double[this.maxExp+1][this.dataSizes.length];
		
		double[][] xaxis = new double[this.maxExp+1][this.dataSizes.length];
		for (int i = 0; i < xaxis.length; i++) {
			for (int j = 0; j < this.dataSizes.length; j++) {
				xaxis[i][j] = this.dataSizes[j];
			}
		}
		
		int baseSize;
		for (int c = 0; c <= maxExp; c++) {
			baseSize = (int) Math.pow(2, c);
			System.out.print("Base: ");
			System.out.print(baseSize);
			System.out.print(", ");
			
			double[][] errors;
			double[] argMinArray = new double[this.dataSizes.length];
			double[] errorMinArray = new double[this.dataSizes.length];
			double truQuery, empQuery, error;
			
			HashMap<String, Matrix>[] empModels;
			HashMap<String, Matrix> emp;
			for (int z = 0; z < repeats; z++){
				errors = new double[this.dataSizes.length][this.maxStates];
				for (int i = 0; i < this.dataSizes.length; i++){
					empModels = this.h.singledataSpectralEmpericalALLMODELS(this.hSize, this.dataSizes[i], this.basisSize, this.maxStates);
					for (int j = 0; j < this.maxStates; j++){		
						//emp = this.h.singledataSpectralEmperical(this.hSize, this.dataSizes[i], this.basisSize, j);
						emp = empModels[j];
						
						for (int q = 0; q < this.maxQuery; q++){
							empQuery = Analysis.probabilityQuery(emp, emp.get("a0"),  emp.get("ainf"), q, baseSize, 2, true);
							truQuery = Analysis.probabilityQuery(this.tru, this.tru.get("a0"),  this.tru.get("ainf"), q, 1, 2, true);
							error = computeError(truQuery, empQuery);
							errors[i][j] += error;
						}
					}
					argMinArray[i] += Analysis.getArgMin( errors[i] ) + 1;
					errorMinArray[i] += Analysis.getMinValue( errors[i] );
				}
				
			}
			
			for (int i = 0; i < errorMinArray.length; i++) {
				argMinArray[i] /= repeats;
				errorMinArray[i] /= repeats; 	
			}
			
			plotErrors[c] = errorMinArray;
			plotArgForErrors[c] = argMinArray;
			
		}
		
		Matrix printBestBaseErrors = new Matrix(plotErrors);
		Matrix printBestBaseArg = new Matrix(plotArgForErrors);
		
		System.out.println("");
		System.out.println("Model Size Errors");
		System.out.println();
		System.out.println("Downwards: BASE, SideWays: #DATA");
		printBestBaseErrors.print(5, 5);
		System.out.println("Downwards: BASE, SideWays: #DATA");
		printBestBaseArg.print(5, 5);
	
		Analysis.outputData(pltFolder + "MinError_Dif_Bases", "X: Data, Y:Min_over_#states", "", xaxis, plotErrors);
		Analysis.outputData(pltFolder + "ArgMin_Dif_Bases", "X: Data, Y:ArgMin_over_#states", "", xaxis, plotArgForErrors);
	}
	
	public void plotTrialsModelSize(int num_lines, int amountOfData){
	
		
		double[][] xaxis = new double[num_lines][this.maxStates];
		double[][] yaxis = new double[num_lines][this.maxStates];
		
		HashMap<String, Matrix>[] empModels;
		HashMap<String, Matrix> emp;
		double truQuery, empQuery, error;
		for (int i = 0; i < num_lines; i++) {
			empModels = this.h.singledataSpectralEmpericalALLMODELS(this.hSize, amountOfData, this.basisSize, xaxis[0].length);
			for (int j = 0; j < empModels.length; j++) {
				emp = empModels[j];
				error = 0;
				for (int c = 0; c < this.maxQuery; c++) {
					truQuery = Analysis.probabilityQuery(this.tru, this.tru.get("a0"),  this.tru.get("ainf"), c, 1, 2, true);
					empQuery = Analysis.probabilityQuery(emp, emp.get("a0"),  emp.get("ainf"), c, 1, 2, true);
					error += computeError(truQuery, empQuery);
				}
				xaxis[i][j] = j+1;
				yaxis[i][j] = error;
			}
			
		}
		
		Analysis.outputData(pltFolder + "Multiple_Trials_ModelError", "X: ModelSize Y:Error", "", xaxis, yaxis);
	}
	
	
	public static void outputData(String filename, String xaxisLabel, String yaxisLabel, double[][] xaxis, double[][] yaxis){
		try {
			PrintWriter writer = new PrintWriter(filename, "UTF-8");
			
			writer.println(xaxisLabel + "," + yaxisLabel);
			
			StringBuilder line = new StringBuilder();
			for (int j = 0; j < xaxis[0].length; j++) {	
				for (int i = 0; i < xaxis.length; i++) {
					line.append(xaxis[i][j] + ",");
					line.append(yaxis[i][j] + " ");
				}
				writer.println( line );
				line.setLength(0);
			} 
			
			writer.close();
			
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	
	public static double getMinValue( double[] a){
		double min = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || min > a[i]){
				init = true;
				min = a[i];
			}
		}
		return min;
	}
	
	public static double getArgMin( double[] a){
		double min = 0;
		double argmin = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || min > a[i]){
				init = true;
				argmin = i;
				min = a[i];
			}
		}
		return argmin;
	}
	
	public static double getMaxValue( double[] a){
		double max = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || max < a[i]){
				init = true;
				max = a[i];
			}
		}
		return max;
	}
	
	public static double getArgMax( double[] a){
		double max = 0;
		double argmax = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || max < a[i]){
				init = true;
				argmax = i;
				max = a[i];
			}
		}
		return argmax;
	}
	
	public static Matrix matrixPower(Matrix m, int exp){
		
		Matrix I = Matrix.identity(m.getArray().length, m.getArray().length);
		if (exp == 0){
			return I;
		}
		else{
			if (exp % 2 == 1) {
				return m.times( matrixPower(m, (exp-1)/2) );
			}
			else{
				return matrixPower(m, exp/2);
			}
		}
	}
	
	public static Matrix matrixQuery(HashMap<String, Matrix> d, int power, int maxPower, int base, boolean forward){	
		int p = maxPower;
	
		int size = d.get( Integer.toString(p) ).getArrayCopy().length;
		Matrix r = Matrix.identity(size,size);
		while(power != 0){
	
			while (p > power){
				p = p/base;
			}
			if (forward){
				r = r.times( d.get(Integer.toString( p ) ) );
			}
			else{
				r = d.get(Integer.toString( p )).times(r);
			}
			power -= p;
		}
		
		return r;
	}
	
	public static double probabilityQuery(HashMap<String, Matrix> d, Matrix ao, Matrix ainf,  int power, int maxPower, int base, boolean forward){
		int p = maxPower;
		Matrix r;
		if (forward){
			r = ao;
		}
		else{
			r = ainf;
		}
		while(power != 0){
			
			while (p > power){
				p = p/base;
			}
			if (forward){
				//System.out.println(p);
				//System.out.println(d.keySet());
				r = r.times( d.get(Integer.toString( p ) ) );
				
			}
			else{
				r = d.get(Integer.toString( p )).times(r);
			}
			power -= p;
		}
		if (forward){
			r = r.times(ainf);
		}
		else{
			r = ao.times(r);
		}
		return r.get(0,0);
		
	}
	
	
	public static Matrix alphaKQuery(HashMap<String, Matrix> d, Matrix ao, int power, int maxPower, int base){	//Always forward for now
		
		int p = maxPower;
		
		Matrix r = ao;
		while(power != 0){
			
			while (p > power){
				p = p/base;
			}
		
			r = r.times( d.get(Integer.toString( p ) ) );
			power -= p;
		}
		
		return r;
		
	}
	
	public static double sumArray(double[] da){
		double s = 0;
		for(double d: da){
			s += d;
		}
		return s;
	}
	
	
	public static double[] incArray(int length){
		double[] output = new double[length];
		for (int i = 0; i < output.length; i++) {
			output[i] = i;
		}
		return output;
	}

}
