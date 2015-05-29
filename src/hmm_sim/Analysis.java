package hmm_sim;

import java.util.ArrayList;
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
		Analysis a = new Analysis(100,10,1000);
		a.compareSigmaError();
		System.out.println("Power Sigma Done");
		a.compareH_Hbar();
		System.out.println("H Hbar Done");
		a.compareHSigmas();
		System.out.println("SigmaError Done");
		a.compareQueryErrors(true);
		System.out.println("Query comp Done");
		a.compareQueryErrors(false);
		System.out.println("Query no comp Done");

	}
	
	public Analysis(int hSize, int basisSize, int repeats){
		this.h = makeHMM();
		this.tru = h.singledataSpectralTrue(hSize, basisSize);
		
		this.hSize = hSize;
		this.basisSize = basisSize;
		
		for (int i = 0; i < repeats; i++) {
			empArray.add(h.singledataSpectralEmperical(hSize, 10000, basisSize));
		}
	}
	
	public void compareH_Hbar(){
		
		HMM h = makeHMM();
		
		int[] sizes = new int[]{100,1000,10000,100000};
		
		HashMap<String, Matrix> emp;
		
		Matrix H;
		Matrix Hbar;
		
		int trailsize = sizes.length;
		int repeats = 10;
		
		double[] dataSize = new double[trailsize];
		double[] error = new double[trailsize];
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
			dataSize[i] = Math.log10(sizes[i]);
			error[i] = avgError;
		}
			
		HelperFunctions.outputData(pltFolder + "True_Hankel_vs_Emperical", "Observation Count","F Norm", dataSize, error );
		
	}
	
	public void compareHSigmas(){
		HashMap<String, Matrix> emp;
		
		int maxQuery = (int) (tru.get("max").norm1()-1); 
		
		double[] sigmaNumber = new double[maxQuery];
		double[] errors = new double[maxQuery];
		
		Matrix h_sigma_true, h_sigma_exp, r;
		int pow;
		
		for (int j = 0; j < empArray.size(); j++) {
			emp = empArray.get(j);
			for (int i = 0; i < errors.length; i++) {
				pow = (int) Math.pow(2, i);
				h_sigma_true = tru.get( Integer.toString(pow) );
				h_sigma_exp = emp.get( Integer.toString(pow) );
				
				r = h_sigma_true.minus( h_sigma_exp );
				
				errors[i] += r.normF();
			}
		}
		
		for (int i = 0; i < maxQuery; i++) {
			sigmaNumber[i] = i;
			errors[i] = errors[i]/empArray.size();
		}
		
		HelperFunctions.outputData(pltFolder + "True_Hx_vs_Emperical_Hx", "Sigma^x","F Norm", sigmaNumber, errors );
		// Add file containing error analysis for alphaInf and alpha0?
	}
	
	public void compareSigmaError(){
		int m = (int) (tru.get("max").norm1()-1); 
		
		double[] sigmaNumber = new double[m];
		double[] errors = new double[m];
		
		Matrix temp1, temp2,  r;
		int pow;
		
		for (int j = 0; j < empArray.size(); j++) {
				
			for (int i = 0; i < m; i++) {
				pow = (int) Math.pow(2, i);
				temp1 = empArray.get(j).get( Integer.toString( pow ) );
				temp1 = HelperFunctions.matrixPower( temp1 , 2);
				temp2 = empArray.get(j).get( Integer.toString(pow*2) );
				r = temp2.minus( temp1 ) ;	
			
				errors[i] += r.normF();
				
			}
			
		}
		
		for (int i = 0; i < errors.length; i++) {
			sigmaNumber[i] = i;
			errors[i] = errors[i]/empArray.size();
		}
		
		HelperFunctions.outputData(pltFolder + "(Hx)^2_v.s H(x^2)", "Sigma^x","F Norm", sigmaNumber, errors );

	}
	
	public void compareQueryErrors(boolean standard){
		int maxexp = (int) tru.get("max").norm1(); //same for tru by construction
		int maxquery = (int) Math.pow(2, maxexp);
		
		double[] querys = new double[maxquery];
		double[] errors = new double[maxquery];
		
		Matrix empQ, truQ;
		Matrix a0emp, ainfemp, a0tru, ainftru, truProb, empProb, error;
		
		HashMap<String, Matrix> emp;
		for (int i = 0; i < maxquery ; i++) {
			truQ = HelperFunctions.matrixQuery(tru, i, 2);
			a0tru = tru.get("a0");
			ainftru = tru.get("ainf").transpose();
			truProb = a0tru.times(truQ).times(ainftru);
			for (int j = 0; j < empArray.size(); j++) {	
				emp = empArray.get(j);
				if(standard){
					empQ = HelperFunctions.matrixQuery(emp, i, 2);
				}
				else{
					empQ = HelperFunctions.matrixPower(emp.get("1"), i);										//inefficient, if slow optimize later
				}
				a0emp = emp.get("a0");
				ainfemp = emp.get("ainf").transpose();
				
				empProb = a0emp.times(empQ).times(ainfemp);
				
				error = truProb.minus(empProb);
				
				errors[i] += Math.abs(error.get(0, 0));
			}
			errors[i] = errors[i]/empArray.size();
			querys[i] = i;
		}
		if (standard){
			HelperFunctions.outputData(pltFolder + "Query_Errors", "Sigma^x","Error", querys, errors );
		}
		else{
			HelperFunctions.outputData(pltFolder + "Query_Errors_Naive", "Sigma^x","Error", querys, errors );
		}
	}
	
	public HMM makeHMM(){
		double[][] p = { {0}, {1}};
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
