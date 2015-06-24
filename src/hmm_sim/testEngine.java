package hmm_sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class testEngine{
	
	public String pltFolder;
	private String fileNameOfDataSet;
	private String[] fileNames;

	private HashMap<Integer, HankelSVDModel[]> empericalModels;
	private QueryEngine[][] fixedSizeQueryEngines;
	private HashMap<Integer,QueryEngine[][]> anySizeQueryEngines;
	private int[] keySetSorted;
	
	private HankelSVDModel trueModel;
	private QueryEngine trueQueryEngine;

	private int base;
	private int basisSize;
	private int maxQuery;
	private int maxStates;
	private int REPEATS;
	
	public static void main(String[] args){
 	
		testEngine a = new testEngine("Models_Emperical_19_12_Toy_Labyrinth/", 40, 2, 50);
	}
	
	public testEngine(String fileNameOfDataSet, int basisSize, int base, int numberPerTrajectorySize){
		this.fileNameOfDataSet = fileNameOfDataSet;
		File[] f = FlowControl.getFiles(fileNameOfDataSet);
		this.fileNames = new String[f.length];
		for (int i = 0; i < fileNames.length; i++) {
			this.fileNames[i] = f[i].getName();
		}
		
		this.basisSize = basisSize;
		this.base = base;
		
		this.trueModel = this.readTrueModel(fileNameOfDataSet);
		this.maxStates = this.trueModel.getRank();
		this.trueQueryEngine = this.trueModel.buildHankelBasedModel(basisSize, base, this.maxStates);
		this.maxQuery = this.trueQueryEngine.getMaxPower() * base;
		
		this.empericalModels = this.readEmpericalModels(fileNameOfDataSet, numberPerTrajectorySize);
				
		this.anySizeQueryEngines = this.getAllSizesQueryEngines(100);
		
		this.keySetSorted = new int[this.anySizeQueryEngines.size()];
		Set<Integer> ks = this.anySizeQueryEngines.keySet();
		int c = 0;
		for (Integer integer : ks) {
			this.keySetSorted[c] = (int) integer;
			c++;
		} 
		Arrays.sort(this.keySetSorted);
		
		for (QueryEngine[][] q: this.anySizeQueryEngines.values()) {
			REPEATS = q.length;
			break;
		}
		
		this.fixedSizePlots();
		System.out.println("Done Fixed Plots");
		
		this.plotBaseDifferences( );
		System.out.println("Done Base Differences");
		
		this.sizeOfModelPlots( );
		System.out.println("Done Model Differences");
		

	}
	
	private HashMap<Integer, QueryEngine[][]> getAllSizesQueryEngines(int numberOfTrajectoriesFromEachSize){
		
		HashMap<Integer, QueryEngine[][]> dataSizeToModels = null;
		QueryEngine[][] Q = new QueryEngine[numberOfTrajectoriesFromEachSize][this.maxStates];
		QueryEngine q;

		try{
			for (String file: this.fileNames) {
				int trajectoryLength = testEngine.getTrajectoryLengthFromFileName(file);
				ObjectInputStream ois = new ObjectInputStream( new FileInputStream(file) );
				int numberOfLinesToRead = numberOfTrajectoriesFromEachSize*3;
				double[] probabilities;
				int basisSize; 
				SingularValueDecomposition svd;
				HankelSVDModel h;
				
				for (int i = 0; i < numberOfLinesToRead; i++) {
					probabilities = (double[]) ois.readObject();
					basisSize = (int) ois.readObject();
					svd = (SingularValueDecomposition) ois.readObject();
					h = new HankelSVDModel(probabilities, basisSize, svd);
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
		HashMap<Integer, HankelSVDModel[]> h = new HashMap<Integer, HankelSVDModel[]>();
		
		HankelSVDModel[] empericalModelsFixedTrajectories;
		try{
			for (String file: fileNames) {	//
				int numberOfTrajectories = testEngine.getTrajectoryLengthFromFileName(file);
				ObjectInputStream ois = new ObjectInputStream(new FileInputStream(file));
				
				empericalModelsFixedTrajectories = new HankelSVDModel[numberPerTrajectorySize];
				for (int i = 0; i < empericalModelsFixedTrajectories.length; i++) {
					empericalModelsFixedTrajectories[i] = (HankelSVDModel) ois.readObject();
				}
				h.put(numberOfTrajectories, empericalModelsFixedTrajectories);
			}
			return h;
		}
		catch(Exception e){
			System.out.println("Problem reading Emperical Models");
			e.printStackTrace();
			return null;
		}
		
	}
	
	public void fixedSizePlots(){		
		
		int maxAhead = 50;
		this.conditionalPlots(maxAhead);
		
		this.compareSigmaError();
		
		this.compareH_Hbar(5);
		
		this.compareASigmas();
		
		System.out.println("SigmaError Done");
		this.compareQueryErrors();
		System.out.println("Query Done");
		
		int num_lines = 10;
		int amountOfData2 = 100;
		this.plotTrialsModelSize(num_lines, amountOfData2);
		System.out.println("Done Multiple Lines Model Error");
	}
	
	public void plotBaseDifferences(){
				
		double[][] dataSize = new double[this.trueQueryEngine.getMaxExponent()][this.anySizeQueryEngines.size()];
		double[][] errors = new double[this.trueQueryEngine.getMaxExponent()][this.anySizeQueryEngines.size()];
		
		double error = 0, empProb, truProb;
		int maxQuery = this.trueQueryEngine.getMaxPower();
		int maxExponent;
		
		for (int c = 0; c < this.keySetSorted.length; c++) {
			QueryEngine[][] qe = this.anySizeQueryEngines.get(keySetSorted[c]);
			for (int i = 0; i < this.REPEATS; i++){
				for (int j = 0; j <= this.trueQueryEngine.getMaxExponent() ; i++){
					maxExponent = (int) Math.pow(this.base, i);
					dataSize[j][c] = this.keySetSorted[c];
					for (int query = 0; query < maxQuery; query++){
						truProb = this.trueQueryEngine.probabilityQuery(query, maxExponent, this.base, true);
						empProb = qe[i][j].probabilityQuery(query, maxExponent, this.base, true);
						error = computeError(truProb, empProb);
						errors[j][c] += error;
					}
				}
			}
			c++;
		}
		
		for (int i = 0; i <= this.trueQueryEngine.getMaxExponent(); i++) {
			for (int j = 0; j < this.anySizeQueryEngines.size(); j++) {
				errors[i][j] /= REPEATS;
			}
		}
		testEngine.outputData(pltFolder + "BaseComp_Area", "X:#Data Seen Y:Fnorm","", dataSize, errors );

		
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

	public void conditionalPlots(QueryEngine[] chosenSizeQueryEngine, int maxK,int maxAhead){
				
		double[][] xaxis = new double[maxK][maxAhead];
		double[][] queryArrayTru = new double[maxK][maxAhead];
		for (int i = 0; i < maxK; i=i++) {
			queryArrayTru[i] = this.trueQueryEngine.conditionalQuery(i, maxAhead);
			xaxis[i] = testEngine.incArray(maxAhead);
		}
		Matrix truPredictions = new Matrix(queryArrayTru);

		double[][] queryArrayEmp = new double[maxK][maxAhead];
		Matrix queryEmpAvg = null, qE = null;
		double[][] errorArray;
		Matrix error, errorAbs, avgError = null;
		
		for (int i = 0; i < this.REPEATS; i++) {
			for (int j = 0; j < maxK; j++) {
				queryArrayEmp[j] = chosenSizeQueryEngine[i].conditionalQuery(j, maxAhead);
			}
			qE = new Matrix(queryArrayEmp);
			if (queryEmpAvg != null){
				queryEmpAvg = queryEmpAvg.plus(qE);
			}
			else{
				queryEmpAvg = qE;
			}
			
			error = qE.minus(truPredictions);
			
			errorArray = new double[maxK][maxAhead];
			for (int j = 0; j < maxK; j++) {
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
		
		avgError = avgError.times(1.0/this.REPEATS);
		queryEmpAvg = queryEmpAvg.times(1.0/this.REPEATS);
		
		testEngine.outputData(pltFolder + "ConditionalError", "x:Traj Length y:|f_k(x)-fhat_k(x)|", "", xaxis, avgError.getArrayCopy());
		testEngine.outputData(pltFolder + "ConditionalEmp", "x:Traj Length y:fhat_k(x)", "", xaxis, queryEmpAvg.getArrayCopy());
		testEngine.outputData(pltFolder + "ConditionalTrue", "x:Traj Length y:f_k(x)", "", xaxis, truPredictions.getArrayCopy());
	}
	
	public void compareH_Hbar(int repeats){
		
		Matrix H = this.trueQueryEngine.getH();
		Matrix Hbar;
		
		double[][] dataSize = new double[1][this.anySizeQueryEngines.size()];
		double[][] error = new double[1][this.anySizeQueryEngines.size()];
		double avgError, e;
	
		int c = 0;
		for (Integer key: this.anySizeQueryEngines.keySet()) {
			QueryEngine[][] q = this.anySizeQueryEngines.get(key);
			avgError = 0;
			for (int j = 0; j < this.REPEATS; j++) {
				Hbar = q[j][0].getH();	//Could replace with HankelSVDMODEL stored globally
				
				e = H.minus(Hbar).normF();
				avgError += e;
			}
			avgError /= this.REPEATS;
			dataSize[0][c] = key;
			error[0][c] = avgError;
			
			c++;
		}
			
		testEngine.outputData(pltFolder + "True_H_vs_Emp", "X:#Data Seen Y:Fnorm","", dataSize, error );
	
	}
	
	public void compareASigmas(QueryEngine[] chosenSizeQueryEngine){
		QueryEngine q;
		
		double[][] sigmaNumber = new double[1][this.trueQueryEngine.getMaxExponent()];
		double[][] errors = new double[1][this.trueQueryEngine.getMaxExponent()];
		
		Matrix a_sigma_true, a_sigma_exp, r;
		
		Matrix[] trueSigmas = this.trueQueryEngine.getAsigmas();
		
		for (int i = 0; i < chosenSizeQueryEngine.length; i++) {
			q = chosenSizeQueryEngine[i];
			Matrix[] empericalSigmas = q.getAsigmas();
			for (int j = 0; j < this.trueQueryEngine.getMaxExponent(); j++) {
				int pow = (int) Math.pow(this.base, j);
				a_sigma_true = trueSigmas[j];
				a_sigma_exp = empericalSigmas[j];

				r = a_sigma_true.minus( a_sigma_exp );
				errors[0][i] += r.norm1();
			}
		}
		
		Matrix[] sigmas = this.trueQueryEngine.getAsigmas();
		for (int i = 0; i < this.trueQueryEngine.getMaxExponent(); i++) {
			int pow = (int) Math.pow(this.base, i);
			sigmaNumber[0][i] = pow;
			errors[0][i] /= (chosenSizeQueryEngine.length*trueSigmas[i].normF());
		}
		
		testEngine.outputData(pltFolder + "True_Ax_vs_Emp", "X:Sigma Y:(T_Ax-E_Ax).Fnorm/T_Ax.norm1","", sigmaNumber, errors );
		// Add file containing error testEngine for alphaInf and alpha0?
	}
	
	public void compareSquareSigmaError(QueryEngine[] chosenSizeQueryEngine){
		
		int maxExpSquareComparison = this.trueQueryEngine.getMaxExponent()-1;
		
		double[][] sigmaNumber = new double[1][maxExpSquareComparison];
		double[][] errors = new double[1][maxExpSquareComparison];
		
		Matrix temp1, temp2, r;
		int pow;
		
		for (int i = 0; i < chosenSizeQueryEngine.length; i++) {
			Matrix[] experimentalSigmas = chosenSizeQueryEngine[i].getAsigmas();
			for (int j = 0; j < maxExpSquareComparison; j++) {
				pow = (int) Math.pow(2, j);
				temp1 = testEngine.matrixPower( experimentalSigmas[j] , this.base);
				temp2 = experimentalSigmas[i+1];
				r = temp2.minus( temp1 );	
				errors[0][j] += r.norm1();
			}
		}
		
		Matrix h_sigma_true;
		Matrix[] trueSigmas = this.trueQueryEngine.getAsigmas();
		for (int j = 0; j < maxExpSquareComparison; j++) {
			pow = (int) Math.pow(2, j+1);
			h_sigma_true = trueSigmas[j];
			sigmaNumber[0][j] = pow;
			errors[0][j] /= (chosenSizeQueryEngine.length*h_sigma_true.norm1());
		}
		
		testEngine.outputData(pltFolder + "(Ax)^2_v.s A(x^2)", "X:Sigma Y:(T_Ax-E_Ax).Fnorm/T_Ax.Fnorm","", sigmaNumber, errors );

	}
	
	public void compareQueryErrors(QueryEngine[] chosenSizeQueryEngine){		
		
		double[][] queries = new double[9][this.maxQuery];
		double[][] errors = new double[9][this.maxQuery];
		
		double[][] baseQueries = new double[this.trueQueryEngine.getMaxExponent()][this.maxQuery];
		double[][] x_base_Queries = new double[this.trueQueryEngine.getMaxExponent()][this.maxQuery];
		
		Matrix a0emp, ainfemp, empQF, empQB;
		Matrix a0tru, ainftru;
		double truProbQF, truProbQB, truProbP , empProbQF, empProbQB, empProbP;
		
		HashMap<String, Matrix> emp;
		for (int i = 0; i < this.maxQuery ; i++) {
	
			a0tru = this.trueQueryEngine.getA0();
			ainftru = this.trueQueryEngine.getAinf();
	
			truProbQF = this.trueQueryEngine.probabilityQuery(i, this.trueQueryEngine.getMaxPower(), this.base, true);
			truProbQB =  this.trueQueryEngine.probabilityQuery(i, this.trueQueryEngine.getMaxPower(), this.base, false);
			truProbP = this.trueQueryEngine.probabilityQuery(i, 1, this.base, true);
			//System.out.println(truProbQF - truProbQB);//Always 0 which makes sense
			//System.out.println(truProbP - truProbQF);
			
			for (int j = 0; j < chosenSizeQueryEngine.length; j++) {	
				QueryEngine q = chosenSizeQueryEngine[j];
		
				empQF = q.matrixQuery( i, this.trueQueryEngine.getMaxPower(), this.base, true);
				empQB = q.matrixQuery( i, this.trueQueryEngine.getMaxPower(), this.base, false);
				
				a0emp = emp.get("a0");
				ainfemp = emp.get("ainf");
				
				empProbQF = q.probabilityQuery( i, this.trueQueryEngine.getMaxPower(), this.base, true);
				empProbQB =  q.probabilityQuery( i, this.trueQueryEngine.getMaxPower(), this.base, false);
				empProbP =  q.probabilityQuery(i, 1, 2, true);
				
				errors[0][i] += Math.abs(truProbQF - empProbQF);	//Tru v.s Base 
				errors[1][i] += truProbQF - empProbQF;

				errors[2][i] += Math.abs(truProbQF - empProbP);		//Tru v.s Naive 
				errors[3][i] += truProbQF - empProbP ;
								
				errors[4][i] += Math.abs(empProbQF - empProbQB );	// Comm Error
				errors[5][i] += empProbQF - empProbQB ;
				
				errors[6][i] += empQF.minus(empQB).normF();			// Matrix Comm Error
				
				double pq;
				for (int k = 0; k < baseQueries.length; k++) {
					pq = Math.abs( q.probabilityQuery( i, (int) Math.pow(2,k), 2, true) - truProbQF);
					baseQueries[k][i] += pq;
				}
			}
			
			for (int j = 0; j < x_base_Queries.length; j++){
				baseQueries[j][i] /= (chosenSizeQueryEngine.length );
				x_base_Queries[j] = testEngine.incArray(this.maxQuery);
			}
			
			for (int c = 0; c < 9; c++) {
				errors[c][i] /= (chosenSizeQueryEngine.length );
				queries[c][i] = i;
			}
		}
		
		testEngine.outputData(pltFolder + "Query_Errors_Base", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,0,2), Arrays.copyOfRange(errors,0,2) );
		testEngine.outputData(pltFolder + "Query_Errors_Naive", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,2,4), Arrays.copyOfRange(errors,2,4) );
		testEngine.outputData(pltFolder + "Comm_Query_Error", "X:Sigma Y:a0(A16A1-A1A16)aI","", Arrays.copyOfRange(queries,4,6), Arrays.copyOfRange(errors,4,6) );
		testEngine.outputData(pltFolder + "Comm_Matrix_Error", "X:Sigma Y:(A16A1-A1A16).Fnorm","", Arrays.copyOfRange(queries,6,7), Arrays.copyOfRange(errors,6,7) );
		testEngine.outputData(pltFolder + "Base_Errors","X:Sigma Y: Error" ,"", x_base_Queries, baseQueries);
		
		double[][] ebase = Arrays.copyOfRange(errors,0, 1);
		double[][] enaive = Arrays.copyOfRange(errors, 2, 3);
		double[][] ejoint = new double[][]{ebase[0], enaive[0]};
		
		double[][] qbase = Arrays.copyOfRange(queries, 2, 3);
		double[][] qnaive = Arrays.copyOfRange(queries, 2, 3);
		double[][] qjoint = new double[][]{qbase[0], qnaive[0]};
		testEngine.outputData(pltFolder + "QError_Base_vs_Naive", "X:Sigma Y:|f(x)-fhat(x)|","",qjoint,ejoint  );
		
		/*	//Print out area under curve between naive and base method
		System.out.println("Highest Max-Base = " + Integer.toString(this.maxPower));
		System.out.println( testEngine.sumArray(errors[0]) );
		System.out.println("Naive Max-Base = 1");
		System.out.println( testEngine.sumArray(errors[2]) );
		*/
		
	}
	
	public double[] conditionalQuery(HashMap<String, Matrix> learned, int k, int maxAhead){
		int maxpow = (int) Math.pow(2,learned.get("max").get(0, 0));
		Matrix alpha_0 = learned.get("a0");
		Matrix alpha_inf = learned.get("ainf");
		//Matrix Ak = testEngine.matrixQuery(learned, k, 2, true);
		Matrix alpha_k = testEngine.alphaKQuery(learned, alpha_0, k, maxpow, 2);//alpha_0.times(Ak);
		
		int nstates = alpha_0.getArray()[0].length;
		
		Matrix mid = Matrix.identity(nstates, nstates).minus(learned.get("1") );
		double normalizer = alpha_k.times( mid.inverse() ).times( alpha_inf ).get(0,0);
		double jointProb;
		
		double[] pA = new double[maxAhead];
		int maxbase = 1;
		for (int i = 0; i < pA.length; i++) {
			jointProb = testEngine.probabilityQuery(learned, alpha_k, alpha_inf, i, maxbase ,2, true);
			pA[i] = jointProb/normalizer;			
		}
		
		return pA;
	}
	
	public void sizeOfModelPlots(int repeats, boolean debug){
		
		double[][] plotErrors = new double[this.trueQueryEngine.getMaxExponent()+1][this.anySizeQueryEngines.size()];
		double[][] plotArgForErrors = new double[this.trueQueryEngine.getMaxExponent()+1][this.anySizeQueryEngines.size()];
		double[][] xaxis = new double[this.trueQueryEngine.getMaxExponent()+1][this.anySizeQueryEngines.size()];
		
		for (int i = 0; i < xaxis.length; i++) {
			for (Integer key :) {
				xaxis[i][c] = key;
				
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
							empQuery = testEngine.probabilityQuery(emp, emp.get("a0"),  emp.get("ainf"), q, baseSize, 2, true);
							truQuery = testEngine.probabilityQuery(this.tru, this.tru.get("a0"),  this.tru.get("ainf"), q, 1, 2, true);
							error = computeError(truQuery, empQuery);
							errors[i][j] += error;
						}
					}
					argMinArray[i] += testEngine.getArgMin( errors[i] ) + 1;
					errorMinArray[i] += testEngine.getMinValue( errors[i] );
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
	
		testEngine.outputData(pltFolder + "MinError_Dif_Bases", "X: Data, Y:Min_over_#states", "", xaxis, plotErrors);
		testEngine.outputData(pltFolder + "ArgMin_Dif_Bases", "X: Data, Y:ArgMin_over_#states", "", xaxis, plotArgForErrors);
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
					truQuery = testEngine.probabilityQuery(this.tru, this.tru.get("a0"),  this.tru.get("ainf"), c, 1, 2, true);
					empQuery = testEngine.probabilityQuery(emp, emp.get("a0"),  emp.get("ainf"), c, 1, 2, true);
					error += computeError(truQuery, empQuery);
				}
				xaxis[i][j] = j+1;
				yaxis[i][j] = error;
			}
			
		}
		
		testEngine.outputData(pltFolder + "Multiple_Trials_ModelError", "X: ModelSize Y:Error", "", xaxis, yaxis);
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
	
	public static int getTrajectoryLengthFromFileName(String filename){
		return Integer.parseInt(filename.split(":")[1]);
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
