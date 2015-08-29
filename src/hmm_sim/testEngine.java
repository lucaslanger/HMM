package hmm_sim;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class testEngine{
	
	public String pltFolder;
	private String[] fileNames;

	private QueryEngine[][] fixedSizeQueryEngines;
	
	HashMap<Integer, QueryEngine[]> fixedModelQE;
	HashMap<Integer, QueryEngine[]> trueRankQueryEngines;
	
	private HashMap<Integer,QueryEngine[][]> anySizeQueryEngines;
	private int[] keySetSorted;
	
	private HankelSVDModel trueModel;
	private HashMap<Integer, HankelSVDModel[]> empericalModels;

	private QueryEngine trueQueryEngine;

	private int base;
	private int basisSize;
	private int maxQuery;
	private int maxStates;
	private int REPEATS;
	private int lowModelSize;
	private int upperModelSize;
	private int numberOfModels;
	private int digitsToPrint;
	private int fixedModelSize;
	private int dataSizeForFixedPlots;
	private int[] modelSizes;
	private HashMap<Integer, Integer> modelSizeToIndex;
	private ModelRetrieval ModelRetrieval;
	private String pltFolderFixed;
	private String pltFolderQError;
	
	public static void main(String[] args){}
	
	public testEngine(String workingFolder, String empModels, String fileNameOfTrueModel, int dataSizeForFixedPlots, int basisSize, int base, int[] modelSizes, int fixedModelSize, int repeats, boolean verbose){
		this.digitsToPrint = 7;
		this.lowModelSize = testEngine.getMinValue(modelSizes);
		this.upperModelSize = testEngine.getMaxValue(modelSizes);
		this.modelSizes = modelSizes;
		this.numberOfModels = modelSizes.length;
		this.fixedModelSize = fixedModelSize;
		
		this.pltFolder = workingFolder + "Plotting_" + empModels + "/";
		this.pltFolderFixed = this.pltFolder + "fixedSize/";
		this.pltFolderQError = this.pltFolder + "qError/";
		testEngine.createFolder(this.pltFolder);
		testEngine.createFolder(this.pltFolderFixed);
		testEngine.createFolder(this.pltFolderQError);
		
		/*this.fileNameOfEmpericalModels = workingFolder + empModels;
		File[] f = testEngine.getFiles(fileNameOfEmpericalModels);
		this.fileNames = new String[f.length];
		for (int i = 0; i < this.fileNames.length; i++) {
			this.fileNames[i] = f[i].getName();
		}
		this.makeKeySetSorted();
		*/

		this.basisSize = basisSize;
		this.base = base;
		this.REPEATS = repeats;
		this.dataSizeForFixedPlots = dataSizeForFixedPlots;
		
		this.ModelRetrieval = new ModelRetrieval(workingFolder, empModels, fileNameOfTrueModel, basisSize, base);
		
		this.keySetSorted = ModelRetrieval.getKeySetSorted();
		this.trueModel = ModelRetrieval.readTrueModel(fileNameOfTrueModel);
		
		int topSingularValuesToPrint = 35;
		System.out.println("Printing out top " + Integer.toString(topSingularValuesToPrint) + " Singular Values!");
		double[] singularValues = putDiagonalToArray(this.trueModel.getSvd().getS(), topSingularValuesToPrint);
		System.out.println( Arrays.toString(singularValues) );
		
		//this.trueModel = this.readTrueModel(workingFolder + fileNameOfTrueModel);
		this.maxStates = this.trueModel.getRank();
		this.modelSizeToIndex = this.initializeModelSizeToIndex();
		this.trueQueryEngine = this.trueModel.buildHankelBasedModel(base, this.maxStates);		
		this.maxQuery = this.trueQueryEngine.getMaxPower()*base; 
	
		double capturedProbability = 0;
	
		for (int i = 0; i <= this.maxQuery; i++) {
			 capturedProbability += this.trueModel.getProbabilities()[i];
		}
		
		if (verbose){
			System.out.println("Hankel Size " + Integer.toString(this.trueModel.getProbabilities().length/2));
			System.out.println("Basis Size: " + Integer.toString(basisSize));
			System.out.println("Repetitions: " + Integer.toString(repeats));
			System.out.println("Base System: " + Integer.toString(base));
			System.out.println("DataSize for Fixed Plots: " + Integer.toString(dataSizeForFixedPlots));
			System.out.println();
	
			System.out.println("True Learned Model Size");
			System.out.println(this.trueModel.getRank());
			System.out.println("Model Range:");
			System.out.print(this.lowModelSize);
			System.out.print("-");
			System.out.println(this.upperModelSize);
	
			System.out.println("Chosen Model = " + Integer.toString(this.fixedModelSize));
			System.out.println("");
			
			System.out.println("MaxQuery:");
			System.out.println(this.maxQuery);
			System.out.println("");
			
			
			System.out.println("Captured Probability: ");
			System.out.println(capturedProbability);
			System.out.println();
			System.out.println("Probabilities:");
			System.out.println( Arrays.toString(this.trueModel.getProbabilities()) );
			System.out.println();
	
		}
		//this.fixedModelQE = ModelRetrieval.getSpecificModelSizeQueryEngines(this.REPEATS, this.fixedModelSize);
		//this.trueRankQueryEngines = ModelRetrieval.getSpecificModelSizeQueryEngines(this.REPEATS, this.trueModel.getRank());
		/*
		if(verbose){
			int topCount = 10;
			double[] e = ModelRetrieval.checkEngine(fixedModelQE.get(dataSizeForFixedPlots)[0], this.trueModel, "FixedModelSize", topCount);
		}*/
		
//		int fixedData = dataSizeForFixedPlots;
		
	}
	
	private double[] putDiagonalToArray(Matrix a, int size){
		double[] r = new double[size];
		for (int i = 0; i < size; i++) {
			r[i] = a.get(i, i);
		}
		return r;
	}
	
	private HashMap<Integer, Integer> initializeModelSizeToIndex(){
		HashMap<Integer, Integer> t = new HashMap<Integer, Integer>();
		for (int i = 0; i < modelSizes.length; i++) {
			t.put(modelSizes[i], i);
		}
		return t;
	}

	static void createFolder(String folder) {
		File dir = new File(folder);
		dir.mkdir();
		
	}

	public void makePlots(){
		
		this.plotBaseDifferences(  );
		System.out.println("Done Base Differences");
		
		ModelEnginePair p = this.ModelRetrieval.getAllSizesQueryEngines(this.REPEATS, this.numberOfModels, this.modelSizes);
		this.anySizeQueryEngines = p.getAnySizeQueryEngines();
		this.empericalModels = p.getEmpericalModels();
		
		this.fixedSizeQueryEngines = this.anySizeQueryEngines.get(this.dataSizeForFixedPlots);
		
		this.fixedSize_Plots();
		System.out.println("Done Fixed Plots");
				
		this.compareH_Hbar();
		System.out.println("Done H, Hbar comparisons");
		
		this.sizeOfModelPlots( );
		System.out.println("Done Model Differences");
		
	}
	
	public void fixedSize_Plots(){		
		
		this.plotTrialsModelSize( this.fixedSizeQueryEngines );

		QueryEngine[] fixedModelSizeEngine = this.fixedModelQE.get( this.dataSizeForFixedPlots );
		
		int maxAhead = this.maxQuery;
		int l = 5;
		int maxbase = 1;
		this.conditionalPlots(fixedModelSizeEngine, l, maxAhead, maxbase);
		
		this.compareSquareSigmaError(fixedModelSizeEngine);
			
		this.compareASigmas();
		
		this.compareQueryErrors(fixedModelSizeEngine);
		
	}
	
	public void plotSingularValues(int numberOfValues){
		Matrix s = this.trueModel.getSvd().getS();
		
		double[][] yaxis = new double[1][numberOfValues];
		double[][] xaxis = new double[1][numberOfValues];
		for (int i = 0; i < numberOfValues; i++) {
			xaxis[0][i] = i;
			yaxis[0][i] = s.get(i,i);
		}
		OutputData.outputData(pltFolder + "SingularValues", "X:ModelSize Y:ith singular value","", xaxis, yaxis);
	}
	
	
	public void modelSizeEffectOverBaseImprovement(String identifier, int fixedDataSize){
		System.out.println("Building Query Engines");
		QueryEngine[][] fixedDataSizeModelEngines = new QueryEngine[this.modelSizes.length][this.REPEATS];
		for (int i = 0; i < fixedDataSizeModelEngines.length; i++) {
			fixedDataSizeModelEngines[i] = this.ModelRetrieval.getSpecificModelSizeQueryEngines(this.REPEATS, this.modelSizes[i]).get(fixedDataSize);
		}
		
		QueryEngine[][] fixedDataCustomBaseEngines = new QueryEngine[this.modelSizes.length][this.REPEATS];
		int numSubstrings = 300;
		int maxBaseSize = 10;
		
		System.out.println("Getting operators");
		System.out.println();
		int[][] operators = this.ModelRetrieval.getOperators(this.REPEATS, maxBaseSize, numSubstrings);
		
		int[] chosenOpForSpeed = operators[0];
		for (int i = 0; i < fixedDataCustomBaseEngines.length; i++) {
			fixedDataCustomBaseEngines[i] = this.ModelRetrieval.getSpecificModelSizeQueryEnginesCustomBase(operators, this.REPEATS, this.modelSizes[i]).get(fixedDataSize);
		}
		
		double[][] errors = new double[this.trueQueryEngine.getMaxExponent()+2][this.modelSizes.length];
		double[][] spreads = new double[this.trueQueryEngine.getMaxExponent()+2][this.modelSizes.length];
		double[][] xAxis = new double[this.trueQueryEngine.getMaxExponent()+2][this.modelSizes.length];
		
		double[] trueP = this.trueModel.getProbabilities();
		
		
		System.out.println("Computing errors!");
		System.out.println();
		for (int r = 0; r < this.REPEATS; r++) {
			for (int j = 0; j < this.trueQueryEngine.getMaxExponent()+1; j++) {
				for (int i = 0; i < this.modelSizes.length; i++) {
					xAxis[j][i] = modelSizes[i];
					double e = 0;
					for (int q = 0; q < this.maxQuery; q++) {
						double p = fixedDataSizeModelEngines[i][r].probabilityQuery(q, (int) Math.pow(2,j), base, true);
						double dif = p - trueP[q];
						e += Math.abs(dif);  
					}
					errors[j][i] += e;
					spreads[j][i] += Math.pow(e, 2);
				}
			}
			
			int lastIndex = this.trueQueryEngine.getMaxExponent()+1;
			for (int i = 0; i < this.modelSizes.length; i++) {
				xAxis[lastIndex][i] = modelSizes[i];
				double e = 0;
				for (int q = 0; q < this.maxQuery; q++) {
					double p = fixedDataCustomBaseEngines[i][r].probabilityQuery(q);
					double dif = p - trueP[q];
					e += Math.abs(dif);  
				}
				errors[lastIndex][i] += e;
				spreads[lastIndex][i] += Math.pow(e, 2);
			}
		}
		
		Matrix ERR = new Matrix(errors).times(1.0/this.REPEATS);
		//ERR.print(5, 5);
		errors = ERR.getArrayCopy();
		
		Matrix ERRSQUARE = new Matrix(spreads).times(1.0/this.REPEATS);
		spreads = ERRSQUARE.getArrayCopy();
		
		for (int j = 0; j < this.trueQueryEngine.getMaxExponent()+2; j++) {
			for (int i = 0; i < this.modelSizes.length; i++) {
				spreads[j][i] = Math.sqrt(spreads[j][i] - Math.pow(errors[j][i],2));
			}
		}
		
		Matrix SPREADS = new Matrix(spreads);
		SPREADS.print(5, 5);
		
		//System.out.println(pltFolder + "BaseImprovementOverModelSizesDatasize:" + Integer.toString(fixedDataSize));
		String title = "Double Loop Timing Predictions";
		String internalComment = "Darker Curves --> Richer Base System";
		System.out.println("Outputting data to: " + identifier);
		
		int L = errors.length;
		int[] rows = new int[]{0,L-2,L-1};
		
		Matrix xT = extractRows(new Matrix(xAxis), rows);
		Matrix yT = extractRows(ERR, rows);
		Matrix sT = extractRows(SPREADS, rows);
		
		yT.print(5, 5);
		System.out.println(pltFolder + identifier);
		OutputData.outputData(pltFolder + identifier, "Model Size", "Error", xT.getArrayCopy(), yT.getArrayCopy(), sT.getArrayCopy(), title, internalComment);
	}
	
	public static Matrix extractRows(Matrix m, int[] rows){
		double[][] r = new double[rows.length][m.getArrayCopy()[0].length];
		int c = 0;
		for (int i : rows) {
			r[c] = m.getArrayCopy()[i];
			c++;
		}
		return new Matrix(r);
	}
	
	
	public void plotBaseDifferences(){			
		double[][] dataSize = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] errors = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] squareErrors = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] variances = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] standardDeviations = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		
		double error = 0, empProb, truProb;
		int maxPower;
		
		for (int c = 0; c < this.keySetSorted.length; c++) {
			QueryEngine[] qe = fixedModelQE.get(this.keySetSorted[c]);
			for (int i = 0; i < this.REPEATS; i++){
				for (int j = 0; j <= this.trueQueryEngine.getMaxExponent() ; j++){
					maxPower = (int) Math.pow(this.base, j);
					dataSize[j][c] = this.keySetSorted[c];
	
					double errorBefore = errors[j][c];
					
					for (int query = 0; query < this.maxQuery; query++){
						truProb = this.trueQueryEngine.probabilityQuery(query, maxPower, this.base, true);
						empProb = qe[i].probabilityQuery(query, maxPower, this.base, true);
						error = this.computeError(truProb, empProb);
						errors[j][c] += error;
					}
					squareErrors[j][c] += Math.pow(errors[j][c] - errorBefore,2);
					
				}
			}
		}
		
		for (int j = 0; j <= this.trueQueryEngine.getMaxExponent(); j++) {
			for (int c = 0; c < keySetSorted.length; c++) {
				errors[j][c] /= this.REPEATS;
				squareErrors[j][c] /= this.REPEATS;
				variances[j][c] = squareErrors[j][c] - Math.pow(errors[j][c],2);
				standardDeviations[j][c] = Math.sqrt(variances[j][c]);
			}
		}
		
		OutputData.outputData(pltFolder + "BaseComp_Area", "X:log(Data) Y:log(Error)","", dataSize, errors );
		
		System.out.println("");
		System.out.println("Base Comp Errors");
		System.out.println("Downwards: BASE, SideWays: DATA");
		Matrix visualErrors = new Matrix(errors);
		visualErrors.print(5, this.digitsToPrint);
		System.out.println("Standard Deviations");
		new Matrix(standardDeviations).print(5, this.digitsToPrint);
		
		System.out.println("Differences");
		new Matrix(testEngine.computeDifferenceWithNaive(errors)).print(5, this.digitsToPrint);
		
		double[][] datasize_differences = testEngine.makeDataSizeDifferences(dataSize);
		
		OutputData.outputData(pltFolder + "Difference Plot_FIXEDMS", "X:log(Data) Y:naive v.s different bases", "", datasize_differences, testEngine.computeDifferenceWithNaive(errors) );

		
	}
	
	private static double[][] makeDataSizeDifferences(double[][] dataSize) {
		double[][] r = new double[dataSize.length-1][dataSize[0].length];
		for (int i = 1; i < dataSize.length; i++) {
			r[i-1] = dataSize[i];
		}
		return r;
	}

	private static double[][] computeDifferenceWithNaive(double[][] d){
		double[][] r = new double[d.length-1][d[0].length];
		double[] t = d[0];
		for (int i = 1; i < d.length; i++) {
			r[i-1] = testEngine.arrayDifference(t, d[i]);
		}
		return r;
	}
	
	private static double[] arrayDifference(double[] t, double[] ds) {
		double[] r = new double[t.length];
		for (int i = 0; i < ds.length; i++) {
			r[i] = t[i] - ds[i];
		}
		return r;
	}

	private double computeError(double truProb, double empProb) {
		double error =  Math.abs( truProb - empProb );
		return error;
	}

	public void conditionalPlots(QueryEngine[] chosenSizeQueryEngine, int maxK, int maxAhead, int maxbase){
				
		double[][] xaxis = new double[maxK][maxAhead];
		double[][] queryArrayTru = new double[maxK][maxAhead];
		for (int i = 0; i < maxK; i++) {
			queryArrayTru[i] = this.conditionalQuery(this.trueQueryEngine ,i, maxAhead, maxbase);
			xaxis[i] = testEngine.incArray(maxAhead);
		}
		Matrix truPredictions = new Matrix(queryArrayTru);

		double[][] queryArrayEmp = new double[maxK][maxAhead];
		Matrix queryEmpAvg = null, qE = null;
		double[][] errorArray;
		Matrix error, errorAbs, avgError = null;
		
		for (int i = 0; i < this.REPEATS; i++) {
			for (int j = 0; j < maxK; j++) {
				queryArrayEmp[j] = this.conditionalQuery(chosenSizeQueryEngine[i], j, maxAhead, maxbase);
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
		
		OutputData.outputData(pltFolder + "ConditionalError", "x:Traj Length y:|f_k(x)-fhat_k(x)|", "", xaxis, avgError.getArrayCopy());
		OutputData.outputData(pltFolder + "ConditionalEmp", "x:Traj Length y:fhat_k(x)", "", xaxis, queryEmpAvg.getArrayCopy());
		OutputData.outputData(pltFolder + "ConditionalTrue", "x:Traj Length y:f_k(x)", "", xaxis, truPredictions.getArrayCopy());
	}
	
	public void compareH_Hbar(){
		
		Matrix H = this.trueModel.getHankel();
		
		double[][] dataSize = new double[1][keySetSorted.length];
		double[][] error = new double[1][keySetSorted.length];
		double avgError, e;
	
		for (int i = 0; i < this.keySetSorted.length; i++) {
			int key = this.keySetSorted[i]; 
			HankelSVDModel[] q = this.empericalModels.get(key);
			
			avgError = 0;
			for (int j = 0; j < this.REPEATS; j++) {
				Matrix Hbar = q[j].getHankel();	
				
				e = H.minus(Hbar).normF();
				avgError += e;
			}
			avgError /= this.REPEATS;
			dataSize[0][i] = Math.log(key);
			error[0][i] = avgError;
			
		}
			
		OutputData.outputData(pltFolder + "True_H_vs_Emp", "X:Data Seen Y:Fnorm","", dataSize, error );
	
	}
			
	
	
	public void compareASigmas(){
		QueryEngine[] chosenSizeQueryEngine = this.trueRankQueryEngines.get(this.dataSizeForFixedPlots);
		
		QueryEngine q;
		
		double[][] sigmaNumber = new double[2][this.trueQueryEngine.getMaxExponent()];
		double[][] errors = new double[2][this.trueQueryEngine.getMaxExponent()];
		
		Matrix a_sigma_true, a_sigma_exp, r;
		
		Matrix[] trueSigmas = this.trueQueryEngine.getAsigmas();
		
		for (int i = 0; i < this.REPEATS; i++) {
			q = chosenSizeQueryEngine[i];
			Matrix[] empericalSigmas = q.getAsigmas();
			for (int j = 0; j < this.trueQueryEngine.getMaxExponent(); j++) {
				a_sigma_true = trueSigmas[j];
				a_sigma_exp = empericalSigmas[j];

				r = a_sigma_true.minus( a_sigma_exp );
				errors[0][j] += r.norm1();
			}
		}
		
		for (int j = 0; j < this.trueQueryEngine.getMaxExponent(); j++) {
			int pow = (int) Math.pow(this.base, j);
			sigmaNumber[0][j] = pow;
			sigmaNumber[1][j] = pow;
			errors[0][j] /= (this.REPEATS*trueSigmas[j].norm1());
			errors[1][j] = (trueSigmas[j].norm1());
		}
		
		OutputData.outputData(pltFolder + "True_Ax_vs_Emp", "X:Sigma Y:(T_Ax-E_Ax).1norm/T_Ax.norm1","", sigmaNumber, errors );
		// Add file containing error testEngine for alphaInf and alpha0?
	}
	
	public void compareSquareSigmaError(QueryEngine[] chosenSizeQueryEngine){
	
		int maxExpSquareComparison = this.trueQueryEngine.getMaxExponent()-1;
		
		double[][] sigmaNumber = new double[3][maxExpSquareComparison];
		double[][] errors = new double[3][maxExpSquareComparison];
		
		Matrix temp1, temp2, temp3, r, r2;
		int pow;
		
		for (int i = 0; i < this.REPEATS; i++) {
			Matrix[] experimentalSigmas = chosenSizeQueryEngine[i].getAsigmas();
			for (int j = 0; j < maxExpSquareComparison; j++) {
				pow = (int) Math.pow(this.base, j);
				temp1 = QueryEngine.matrixPower( experimentalSigmas[j] , this.base);
				temp2 = experimentalSigmas[j+1];
				temp3 = QueryEngine.matrixPower( experimentalSigmas[0], pow*this.base );

				r = temp2.minus( temp1 );	
				r2 = temp2.minus( temp3 );
				errors[0][j] += r.norm1();
				errors[2][j] += r2.norm1();
			}
		}
		
		Matrix h_sigma_true;
		Matrix[] trueSigmas = this.trueQueryEngine.getAsigmas();
		for (int j = 0; j < maxExpSquareComparison; j++) {
			pow = (int) Math.pow(2, j+1);
			h_sigma_true = trueSigmas[j+1];
			sigmaNumber[0][j] = pow;
			sigmaNumber[1][j] = pow;
			sigmaNumber[2][j] = pow;
			errors[0][j] /= (this.REPEATS*h_sigma_true.norm1());
			errors[2][j] /= (this.REPEATS*h_sigma_true.norm1());
			errors[1][j] = h_sigma_true.norm1();
			
		}
		
		OutputData.outputData(pltFolder + "(Ax)^2_v.s A(x^2)", "X:Sigma Y:(T_Ax-E_Ax).1norm/T_Ax.1norm","", sigmaNumber, errors );

	}
	
	public void compareQueryErrors(QueryEngine[] chosenSizeQueryEngine){		
		
		double[][] queries = new double[9][this.maxQuery];
		double[][] errors = new double[9][this.maxQuery];
		
		double[][] baseQueries = new double[this.trueQueryEngine.getMaxExponent()+1][this.maxQuery];
		double[][] x_base_Queries = new double[this.trueQueryEngine.getMaxExponent()+1][this.maxQuery];
		
		double truProbQF,  empProbQF, empProbQB, empProbP;
		//double truProbQB, truProbP;
		//Matrix empQF, empQB;
		for (int i = 0; i < this.maxQuery ; i++) {
	
			truProbQF = this.trueQueryEngine.probabilityQuery(i, this.trueQueryEngine.getMaxPower(), this.base, true);
			//truProbQB =  this.trueQueryEngine.probabilityQuery(i, this.trueQueryEngine.getMaxPower(), this.base, false);
			//truProbP = this.trueQueryEngine.probabilityQuery(i, 1, this.base, true);
			//System.out.println(truProbQF - truProbQB);//Always 0 which makes sense
			//System.out.println(truProbP - truProbQF);
			
			for (int j = 0; j < this.REPEATS; j++) {	
				QueryEngine q = chosenSizeQueryEngine[j];
		
				//empQF = q.matrixQuery( i, this.trueQueryEngine.getMaxPower(), this.base, true);
				//empQB = q.matrixQuery( i, this.trueQueryEngine.getMaxPower(), this.base, false);
				//errors[6][i] += empQF.minus(empQB).normF();			// Matrix Comm Error
				
				empProbQF = q.probabilityQuery( i, this.trueQueryEngine.getMaxPower(), this.base, true);
				empProbQB = q.probabilityQuery( i, this.trueQueryEngine.getMaxPower(), this.base, false);
				empProbP = q.probabilityQuery(i, 1, this.base, true);
				
				errors[0][i] += Math.abs(truProbQF - empProbQF);	//Tru v.s Base 
				errors[1][i] += truProbQF - empProbQF;

				errors[2][i] += Math.abs(truProbQF - empProbP);		//Tru v.s Naive 
				errors[3][i] += truProbQF - empProbP ;
								
				errors[4][i] += Math.abs(empProbQF - empProbQB );	// Comm Error
				errors[5][i] += empProbQF - empProbQB ;
								
				double pq;
				for (int k = 0; k < baseQueries.length; k++) {
					pq = Math.abs( q.probabilityQuery( i, (int) Math.pow(2,k), 2, true) - truProbQF);
					baseQueries[k][i] += pq;
				}
			}
			
			for (int j = 0; j < x_base_Queries.length; j++){
				baseQueries[j][i] /= this.REPEATS;
				x_base_Queries[j] = testEngine.incArray(this.maxQuery);
			}
			
			for (int c = 0; c < 9; c++) {
				errors[c][i] /= this.REPEATS;
				queries[c][i] = i;
			}
		}
		
		/*System.out.println("Commutative Avg error");
		System.out.println( Arrays.toString(errors[4]) );
		System.out.println();
		*/
		
		OutputData.outputData(pltFolder + "Query_Errors_Base", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,0,2), Arrays.copyOfRange(errors,0,2) );
		OutputData.outputData(pltFolder + "Query_Errors_Naive", "X:Sigma Y: Green:Absolute","", Arrays.copyOfRange(queries,2,4), Arrays.copyOfRange(errors,2,4) );
		OutputData.outputData(pltFolder + "Comm_Qerror", "X:Sigma Y:a0(A16A1-A1A16)aI","", Arrays.copyOfRange(queries,4,6), Arrays.copyOfRange(errors,4,6) );
		OutputData.outputData(pltFolder + "Comm_Merror", "X:Sigma Y:(A16A1-A1A16).Fnorm","", Arrays.copyOfRange(queries,6,7), Arrays.copyOfRange(errors,6,7) );
		OutputData.outputData(pltFolder + "Base_Errors","X:Sigma Y: Error" ,"", x_base_Queries, baseQueries);
		
		double[][] ebase = Arrays.copyOfRange(errors, 0, 1);
		double[][] enaive = Arrays.copyOfRange(errors, 2, 3);
		double[][] ejoint = new double[][]{ebase[0], enaive[0]};
		
		double[][] qbase = Arrays.copyOfRange(queries, 2, 3);
		double[][] qnaive = Arrays.copyOfRange(queries, 2, 3);
		double[][] qjoint = new double[][]{qbase[0], qnaive[0]};
		OutputData.outputData(pltFolder + "QError_Base_vs_Naive", "X:Sigma Y:|f(x)-fhat(x)|","", qjoint,ejoint  );
		
	}
	
	public double[] conditionalQuery(QueryEngine q, int k, int maxAhead, int maxbase ){
		int maxpow = (int) Math.pow(this.base,q.getMaxExponent());
		Matrix alpha_k = q.alphaKQuery( k, maxpow, 2);
		
		int nstates = q.getA0().getArray()[0].length;
		
		Matrix mid = Matrix.identity(nstates, nstates).minus(q.getAsigmas()[0] );
		double normalizer = alpha_k.times( mid.inverse() ).times( q.getAinf() ).get(0,0);
		double jointProb;
		
		double[] pA = new double[maxAhead];
	
		for (int i = 0; i < pA.length; i++) {
			jointProb = q.probabilityQuery( i+k, maxbase ,2, true);
			pA[i] = jointProb/normalizer;			
		}
		
		return pA;
	}
	
	public void sizeOfModelPlots(){
		
		double[][] plotErrors = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] plotArgForErrors = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] xaxis = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		
		double[][] squareErrorMin = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] variances = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] standardDeviations = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		
		for (int i = 0; i <= this.trueQueryEngine.getMaxExponent(); i++) {
			for (int j = 0; j < this.keySetSorted.length; j++) {
				xaxis[i][j] = this.keySetSorted[j];
			}
		}
		
		int baseSize;
		for (int c = 0; c <= this.trueQueryEngine.getMaxExponent(); c++) {
			baseSize = (int) Math.pow(this.base, c);
			
			double[][] errors;
			double[] argMinArray = new double[this.keySetSorted.length];
			double[] errorMinArray = new double[this.keySetSorted.length];
			double truQuery, empQuery, error, minValue;
			
			errors = new double[this.keySetSorted.length][this.numberOfModels];
			for (int i = 0; i < this.keySetSorted.length; i++){
				for (int z = 0; z < this.REPEATS; z++){		
					errors[i] = new double[this.numberOfModels];
					for (int j = 0; j < this.numberOfModels; j++){			
						QueryEngine q = this.anySizeQueryEngines.get(this.keySetSorted[i])[j][z];
						for (int k = 0; k < this.maxQuery; k++){
							empQuery = q.probabilityQuery(k, baseSize, this.base, true);
							truQuery = this.trueQueryEngine.probabilityQuery(k, 1, this.base, true);
							error = computeError(truQuery, empQuery);
							errors[i][j] += error;
						}
					}
					argMinArray[i] += modelSizes[ (int) testEngine.getArgMin( errors[i] )];
					
					minValue = testEngine.getMinValue( errors[i] );
					errorMinArray[i] += minValue;
					squareErrorMin[c][i] += Math.pow(minValue,2);
				}
			}
			
			for (int i = 0; i < errorMinArray.length; i++) {
				argMinArray[i] /= this.REPEATS;
				
				errorMinArray[i] /= this.REPEATS;
				squareErrorMin[c][i] /= this.REPEATS;
				variances[c][i] = squareErrorMin[c][i] - Math.pow(errorMinArray[i],2);
				standardDeviations[c][i] = Math.sqrt(variances[c][i]);
			}
			
			plotErrors[c] = errorMinArray;
			plotArgForErrors[c] = argMinArray;
			
		}
		
		Matrix printBestBaseErrors = new Matrix(plotErrors);
		Matrix printBestBaseArg = new Matrix(plotArgForErrors);
		
		System.out.println("");
		System.out.println("Model Size Errors");
		System.out.println();
		System.out.println("Downwards: BASE, SideWays: DATA");
		printBestBaseErrors.print(5, this.digitsToPrint);
		System.out.println("Averaged Best Model Sizes");
		printBestBaseArg.print(5, this.digitsToPrint);
		System.out.println("Standard Deviations");
		new Matrix(standardDeviations).print(5, this.digitsToPrint);
	
		OutputData.outputData(pltFolder + "MinError_Dif_Bases", "X: log(Data) Y:log(Min_over_states)", "", xaxis, plotErrors);
		OutputData.outputData(pltFolder + "ArgMin_Dif_Bases", "X: log(Data) Y:ArgMin_over_states", "", xaxis, plotArgForErrors);
		
		System.out.println("Differences:");
		new Matrix( testEngine.computeDifferenceWithNaive(plotErrors) ).print(5, this.digitsToPrint);
		
		double[][] datasize_differences = testEngine.makeDataSizeDifferences(xaxis);
		
		OutputData.outputData(pltFolder + "Difference Plot", "X:log(Data) Y:naive v.s different bases", "", datasize_differences, testEngine.computeDifferenceWithNaive(plotErrors) );

	}
	
	public void plotTrialsModelSize(QueryEngine[][] chosenSizeQueryEngine ){
		
		double[][] xaxis = new double[this.numberOfModels][this.REPEATS];
		double[][] yaxis = new double[this.numberOfModels][this.REPEATS];
		
		double truQuery, empQuery, error;
		
		for (int i = 0; i < this.numberOfModels; i++) {
			for (int j = 0; j < this.REPEATS; j++) {
				QueryEngine q = chosenSizeQueryEngine[i][j];
				error = 0;
				for (int c = 0; c < this.maxQuery; c++) {
					truQuery = this.trueQueryEngine.probabilityQuery( c, 1, 2, true);
					empQuery = q.probabilityQuery( c, 1, 2, true);
					error += computeError(truQuery, empQuery);
				}
				xaxis[i][j] = this.modelSizes[i];
				yaxis[i][j] = error;
			}
			
		}
		xaxis = new Matrix(xaxis).transpose().getArrayCopy();
		yaxis = new Matrix(yaxis).transpose().getArrayCopy();
		OutputData.outputData(pltFolder + "Multiple_Trials_ModelError", "X: ModelSize Y:Error", "", xaxis, yaxis);
		//System.out.println("Fixed Size repeated trial plots");
		//new Matrix(yaxis).print(5, this.digitsToPrint);
	}
	
	

	
	private static int getMinValue( int[] a){
		int min = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || min > a[i]){
				init = true;
				min = a[i];
			}
		}
		return min;
	}
		
	
	private static int getMaxValue( int[] a){
		int max = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (!init || max < a[i]){
				init = true;
				max = a[i];
			}
		}
		return max;
	}
	
	private static double getMinValue( double[] a){
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
	
	private static double getArgMin( double[] a){
		double min = 0;
		double argmin = 0;
		boolean init = false;
		for (int i = 0; i < a.length; i++) {
			if (init == false || min > a[i]){
				init = true;
				argmin = i;
				min = a[i];
			}
		}
		
		return argmin;
	}
	
	private static double getMaxValue( double[] a){
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
	
	private static double getArgMax(double[] a){
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
	
	private static double sumArray(double[] da){
		double s = 0;
		for(double d: da){
			s += d;
		}
		return s;
	}
	
	
	private static double[] incArray(int length){
		double[] output = new double[length];
		for (int i = 0; i < output.length; i++) {
			output[i] = i;
		}
		return output;
	}
	
	private static int getTrajectoryLengthFromFileName(String filename){
		return Integer.parseInt(filename.split(":")[1]);
	}
	
	
	public static double[][] getTopErrorIndices(double[] r, int topcount){
		double[] d = new double[r.length];
		for (int i = 0; i <r.length; i++) {
			d[i] = r[i];
		}
		HashMap<Double, ArrayList<Integer>> val_to_index = new HashMap<Double, ArrayList<Integer>>();
		for (int i = 0; i < d.length; i++) {
			ArrayList<Integer> v = val_to_index.get(d[i]);
			if (v == null){
				v = new ArrayList<Integer>();
			}
			v.add(i);
			val_to_index.put(d[i], v);
		}
		Arrays.sort(d);
		double[] reverse_d = new double[d.length];
		for (int i = 0; i < d.length; i++) {
			reverse_d[d.length-1-i] = d[i];
		}
		d = reverse_d;
	
		double[] top = new double[topcount];
		double[] topvalues = new double[topcount];
		
		int i = 0;
		int t = 0;
		while(t < top.length){
			while(t-i <val_to_index.get(d[i]).size() && t < top.length){
				top[t] = val_to_index.get(d[i]).get(t-i);
				topvalues[t] = d[i];
				t++;
			}
			i++;
		}
		return new double[][]{ top, topvalues };
		
		
	}
	

}
