package hmm_sim;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Set;

import Jama.Matrix;

public class testEngine{
	
	public String pltFolder;
	private String fileNameOfEmpericalModels;
	private String[] fileNames;

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
	private int lowModelSize;
	private int upperModelSize;
	private int numberOfModels;
	
	public static void main(String[] args){}
	
	public testEngine(String workingFolder, String empModels, String fileNameOfTrueModel, int dataSizeForFixedPlots, int basisSize, int base, int numberPerTrajectorySize){
		this.fileNameOfEmpericalModels = workingFolder + empModels;
		this.pltFolder = workingFolder + "Plotting_" + empModels + "/";
		testEngine.createFolder(this.pltFolder);
		File[] f = testEngine.getFiles(fileNameOfEmpericalModels);
		this.fileNames = new String[f.length];
		for (int i = 0; i < this.fileNames.length; i++) {
			this.fileNames[i] = f[i].getName();
		}
		
		this.basisSize = basisSize;
		this.base = base;
		this.trueModel = this.readTrueModel(fileNameOfTrueModel);
		this.maxStates = this.trueModel.getRank();
		this.lowModelSize = this.trueModel.getRank()-2;
		this.upperModelSize = this.trueModel.getRank()+2;
		this.numberOfModels = this.upperModelSize - this.lowModelSize;
		System.out.println("True Learned Model Size");
		System.out.println(this.trueModel.getRank());
		System.out.println("Model Range:");
		System.out.print(this.lowModelSize);
		System.out.print("-");
		System.out.println(this.upperModelSize);
		
		this.trueQueryEngine = this.trueModel.buildHankelBasedModel(this.basisSize, base, this.maxStates);
		this.maxQuery = this.trueQueryEngine.getMaxPower(); //* base;
		System.out.println("MaxQuery:");
		System.out.println(this.maxQuery);
		this.REPEATS = numberPerTrajectorySize;
		
		this.anySizeQueryEngines = this.getAllSizesQueryEngines(numberPerTrajectorySize);
		this.makeKeySetSorted();
				
		this.fixedSizeQueryEngines = this.anySizeQueryEngines.get(dataSizeForFixedPlots);
		
		this.makePlots();

	}
	

	private static File[] getFiles(String folder) {
		File dir = new File(folder);
		return dir.listFiles();
	}

	private static void createFolder(String folder) {
		File dir = new File(folder);
		dir.mkdir();
		
	}

	private void makePlots(){
		//this.fixedSize_Plots();
		//System.out.println("Done Fixed Plots");
		
		//this.compareH_Hbar(5);
		//System.out.println("Done H, Hbar comparisons");
		
		this.plotBaseDifferences( );
		System.out.println("Done Base Differences");
		
		this.sizeOfModelPlots( );
		System.out.println("Done Model Differences");
		
	}
	
	private void makeKeySetSorted(){
		this.keySetSorted = new int[this.anySizeQueryEngines.size()];
		Set<Integer> ks = this.anySizeQueryEngines.keySet();
		int c = 0;
		for (Integer integer : ks) {
			this.keySetSorted[c] = (int) integer;
			c++;
		} 
		Arrays.sort(this.keySetSorted);
	}
	
	private HashMap<Integer, QueryEngine[][]> getAllSizesQueryEngines(int numberOfTrajectoriesFromEachSize){
		
		HashMap<Integer, QueryEngine[][]> dataSizeToModels = new HashMap<Integer, QueryEngine[][]>();
		QueryEngine[][] Q = new QueryEngine[this.numberOfModels][numberOfTrajectoriesFromEachSize];;
		QueryEngine q;

		try{
			for (String f: this.fileNames) {
				Q = new QueryEngine[this.numberOfModels][numberOfTrajectoriesFromEachSize];	//Weird bug
				String file = this.fileNameOfEmpericalModels + f;
				int trajectoryLength = testEngine.getTrajectoryLengthFromFileName(file);
				ObjectInputStream ois = new ObjectInputStream( new FileInputStream(file) );
				HankelSVDModel h;
				
				for (int i = 0; i < numberOfTrajectoriesFromEachSize; i++) {
					h = (HankelSVDModel) ois.readObject();
					for (int j = 0; j < this.numberOfModels; j++) {
						q = h.buildHankelBasedModel(this.basisSize, this.base, this.lowModelSize + j);
						Q[j][i] = q;
					}
				}
				dataSizeToModels.put(trajectoryLength, Q);
				ois.close();
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
			ois.close();
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
	
	public void fixedSize_Plots(){		
		
		this.plotTrialsModelSize( this.fixedSizeQueryEngines );
		
		int modelSize = this.trueQueryEngine.getAsigmas()[0].getArrayCopy().length;	//Rank of True Model
		QueryEngine[] fixedModelSizeEngine = this.fixedSizeQueryEngines[modelSize-1];
		
		int maxAhead = 50, l = 5;
		this.conditionalPlots(fixedModelSizeEngine, l, maxAhead);
		
		this.compareSquareSigmaError(fixedModelSizeEngine);
			
		this.compareASigmas(fixedModelSizeEngine);
		
		this.compareQueryErrors(fixedModelSizeEngine);
	}
	
	public void plotBaseDifferences(){
		int modelSize = this.numberOfModels/2;		
		
		double[][] dataSize = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] errors = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] squareErrors = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] variances = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];
		double[][] standardDeviations = new double[this.trueQueryEngine.getMaxExponent()+1][this.keySetSorted.length];

		
		double error = 0, empProb, truProb;
		int maxPower;
		
		for (int c = 0; c < this.keySetSorted.length; c++) {
			QueryEngine[] qe = this.anySizeQueryEngines.get(this.keySetSorted[c])[modelSize];
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
		
		testEngine.outputData(pltFolder + "BaseComp_Area", "X:log(Data) Y:log(Error)","", dataSize, errors );
		
		System.out.println("");
		System.out.println("Base Comp Errors Modelsize=" + Integer.toString(this.maxStates));
		System.out.println("Downwards: BASE, SideWays: #DATA");
		Matrix visualErrors = new Matrix(errors);
		visualErrors.print(5, 15);
		System.out.println("Standard Deviations");
		new Matrix(standardDeviations).print(5, 15);
		
	}
	
	private double computeError(double truProb, double empProb) {
		double error =  Math.abs( truProb - empProb );
		return error;
	}

	public void conditionalPlots(QueryEngine[] chosenSizeQueryEngine, int maxK, int maxAhead){
				
		double[][] xaxis = new double[maxK][maxAhead];
		double[][] queryArrayTru = new double[maxK][maxAhead];
		for (int i = 0; i < maxK; i++) {
			queryArrayTru[i] = this.conditionalQuery(this.trueQueryEngine ,i, maxAhead);
			xaxis[i] = testEngine.incArray(maxAhead);
		}
		Matrix truPredictions = new Matrix(queryArrayTru);

		double[][] queryArrayEmp = new double[maxK][maxAhead];
		Matrix queryEmpAvg = null, qE = null;
		double[][] errorArray;
		Matrix error, errorAbs, avgError = null;
		
		for (int i = 0; i < this.REPEATS; i++) {
			for (int j = 0; j < maxK; j++) {
				queryArrayEmp[j] = this.conditionalQuery(chosenSizeQueryEngine[i], j, maxAhead);
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
		
		double[][] dataSize = new double[1][keySetSorted.length];
		double[][] error = new double[1][keySetSorted.length];
		double avgError, e;
	
		int c = 0;
		for (int i = 0; i < this.keySetSorted.length; i++) {
			int key = this.keySetSorted[i]; 
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
				a_sigma_true = trueSigmas[j];
				a_sigma_exp = empericalSigmas[j];

				r = a_sigma_true.minus( a_sigma_exp );
				errors[0][j] += r.norm1();
			}
		}
		
		for (int j = 0; j < this.trueQueryEngine.getMaxExponent(); j++) {
			int pow = (int) Math.pow(this.base, j);
			sigmaNumber[0][j] = pow;
			errors[0][j] /= (chosenSizeQueryEngine.length*trueSigmas[j].normF());
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
				temp1 = QueryEngine.matrixPower( experimentalSigmas[j] , this.base);
				temp2 = experimentalSigmas[j+1];
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
		
		Matrix empQF, empQB;
		double truProbQF,  empProbQF, empProbQB, empProbP;
		//double truProbQB, truProbP;

		for (int i = 0; i < this.maxQuery ; i++) {
	
			truProbQF = this.trueQueryEngine.probabilityQuery(i, this.trueQueryEngine.getMaxPower(), this.base, true);
			//truProbQB =  this.trueQueryEngine.probabilityQuery(i, this.trueQueryEngine.getMaxPower(), this.base, false);
			//truProbP = this.trueQueryEngine.probabilityQuery(i, 1, this.base, true);
			//System.out.println(truProbQF - truProbQB);//Always 0 which makes sense
			//System.out.println(truProbP - truProbQF);
			
			for (int j = 0; j < chosenSizeQueryEngine.length; j++) {	
				QueryEngine q = chosenSizeQueryEngine[j];
		
				empQF = q.matrixQuery( i, this.trueQueryEngine.getMaxPower(), this.base, true);
				empQB = q.matrixQuery( i, this.trueQueryEngine.getMaxPower(), this.base, false);
				
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
	
	public double[] conditionalQuery(QueryEngine q, int k, int maxAhead){
		int maxpow = (int) Math.pow(this.base,q.getMaxExponent());
		Matrix alpha_k = q.alphaKQuery( k, maxpow, 2);
		
		int nstates = q.getA0().getArray()[0].length;
		
		Matrix mid = Matrix.identity(nstates, nstates).minus(q.getAsigmas()[0] );
		double normalizer = alpha_k.times( mid.inverse() ).times( q.getAinf() ).get(0,0);
		double jointProb;
		
		double[] pA = new double[maxAhead];
		int maxbase = 1;
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
			baseSize = (int) Math.pow(2, c);
			
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
					argMinArray[i] += testEngine.getArgMin( errors[i] ) + this.lowModelSize;
					
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
		System.out.println("Downwards: BASE, SideWays: #DATA");
		printBestBaseErrors.print(5, 5);
		System.out.println("Averaged Best Model Sizes");
		printBestBaseArg.print(5, 5);
		System.out.println("Standard Deviations");
		new Matrix(standardDeviations).print(5, 10);
	
		testEngine.outputData(pltFolder + "MinError_Dif_Bases", "X: log(Data) Y:log(Min_over_#states)", "", xaxis, plotErrors);
		testEngine.outputData(pltFolder + "ArgMin_Dif_Bases", "X: log(Data) Y:ArgMin_over_#states", "", xaxis, plotArgForErrors);
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
				xaxis[i][j] = i+this.lowModelSize;
				yaxis[i][j] = error;
			}
			
		}
		xaxis = new Matrix(xaxis).transpose().getArrayCopy();
		yaxis = new Matrix(yaxis).transpose().getArrayCopy();
		testEngine.outputData(pltFolder + "Multiple_Trials_ModelError", "X: ModelSize Y:Error", "", xaxis, yaxis);
		System.out.println("Fixed Size repeated trial plots");
		new Matrix(yaxis).print(5, 5);
	}
	
	
	private static void outputData(String filename, String xaxisLabel, String yaxisLabel, double[][] xaxis, double[][] yaxis){
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
			e.printStackTrace();
		}

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
	
	private static double getArgMax( double[] a){
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

}
