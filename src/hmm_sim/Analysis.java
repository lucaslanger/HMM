package hmm_sim;

import java.util.HashMap;

import Jama.Matrix;

public class Analysis {
	
	public static final String pltFolder = "plotting/";
	
	public static void main(String[] args){
		Analysis a = new Analysis();
		//a.compareSigmaError();
		//a.compareH_Hbar();
		a.compareHSigmas();
	}
	
	public Analysis(){
		
	}
	
	public void compareH_Hbar(){
		
		HMM h = makeHMM();
		
		int[] sizes = new int[]{100,1000,10000,100000};
		int repeats = 10;
		int basisSize = 5;
		int hSize = 100;
		
		HashMap<String, Matrix> emp;
		HashMap<String, Matrix> tru = h.singledataSpectralTrue(hSize, basisSize);
		
		Matrix H;
		Matrix Hbar;
		
		int trailsize = sizes.length;
		double[] dataSize = new double[trailsize];
		double[] error = new double[trailsize];
		double avgError, e;
		
		for (int i = 0; i < trailsize; i++) {
			//System.out.println("Trial size");
			//System.out.println(sizes[i]);
			avgError = 0;
			for (int j = 0; j < repeats; j++) {
				emp = h.singledataSpectralEmperical(hSize, sizes[i], basisSize);
				
				H = tru.get("H");
				Hbar = emp.get("H");
				
				e = H.minus(Hbar).normF();
				avgError += e;
			}
			avgError /= repeats;
			dataSize[i] = sizes[i];
			error[i] = avgError;
		}
			
		HelperFunctions.outputData(pltFolder + "H_error", "Sample size","Error", dataSize, error );
		
	}
	
	public void compareHSigmas(){
		HMM h = makeHMM();
		HashMap<String, Matrix> emp = h.singledataSpectralEmperical(100,100000,10);
		HashMap<String, Matrix> tru = h.singledataSpectralTrue(100,10);
		
		compareQueryErrors(emp, tru);
		/*
		System.out.println("H error");
		emp.get("H").minus(tru.get("H")).print(5,5);
		
		System.out.println("P's");
		emp.get("pinv").print(5, 5);
		tru.get("pinv").print(5, 5);
		emp.get("pinv").minus(tru.get("pinv")).print(5, 5);
		
		System.out.println("SVDs");
		emp.get("s_values").print(5, 5);
		tru.get("s_values").print(5, 5);
		
		emp.get("U").print(5, 5);
		tru.get("U").print(5, 5);
		
		emp.get("VT").print(5, 5);
		tru.get("VT").print(5, 5);
		
		System.out.println("Hsigma=1 error");
		emp.get("1").minus(tru.get("1")).print(5,5);
		*/
		
		int m = (int) (emp.get("max").norm1()-1); 
		
		double[] sigmaNumber = new double[m];
		double[] error = new double[m];
		
		Matrix h_sigma_true, h_sigma_exp, r;
		int pow;
		for (int i = 0; i < error.length; i++) {
			pow = (int) Math.pow(2, i);
			h_sigma_true = tru.get( Integer.toString(pow) );
			h_sigma_exp = emp.get( Integer.toString(pow) );
			
			//h_sigma_exp.print(5, 5);
			//h_sigma_true.print(5, 5);
			
			r = h_sigma_true.minus( h_sigma_exp );
			
			sigmaNumber[i] = pow;
			error[i] = r.normF();
		}
		
		HelperFunctions.outputData(pltFolder + "Sigma_Htrue_Hbar", "Sigma","Error", sigmaNumber, error );

	}
	
	public static void compareQueryErrors(HashMap<String, Matrix> emp, HashMap<String, Matrix> tru){
		int maxexp = (int) emp.get("max").norm1(); //same for tru by construction
		int maxquery = (int) Math.pow(2, maxexp);
		
		double[] querys = new double[maxquery];
		double[] errors = new double[maxquery];
		
		Matrix empQ, truQ;
		Matrix a0emp, ainfemp, a0tru, ainftru, error;
		for (int i = 0; i < maxquery ; i++) {
			empQ = HelperFunctions.matrixQuery(emp, i, 2);
			truQ = HelperFunctions.matrixQuery(tru, i, 2);
			
			//System.out.println("Error between emp and tru for sigma^" + Integer.toString(i));
			a0emp = emp.get("a0");
			ainfemp = emp.get("ainf").transpose();
			
			a0tru = tru.get("a0");
			ainftru = tru.get("ainf").transpose();
			
			error = (a0tru.times(truQ).times(ainftru)).minus( a0emp.times(empQ).times(ainfemp) );
			querys[i] = i;
			errors[i] = error.norm1();
		}
		
		HelperFunctions.outputData(pltFolder + "Query_Errors", "Sigma","Error", querys, errors );
		
	}
	
	public void compareSigmaError(){
		HMM h = makeHMM();
		HashMap<String, Matrix> emp = h.singledataSpectralEmperical(100,100000,10);
		int m = (int) (emp.get("max").norm1()-1); 
		
		double[] sigmaNumber = new double[m];
		double[] error = new double[m];
		
		Matrix temp1, temp2,  r;
		int pow;
		for (int i = 0; i < m; i++) {
			pow = (int) Math.pow(2, i);
			temp1 = emp.get( Integer.toString( pow ) );
			temp1 = HelperFunctions.matrixPower( temp1 , 2);
			temp2 = emp.get( Integer.toString(pow*2) );
			r = temp2.minus( temp1 ) ;	
		
			sigmaNumber[i] = pow;
			error[i] = r.normF();
			/*
			System.out.println("Error between consecutive Asigmas");
			System.out.println(pow);
			System.out.println(r.normF());
			System.out.println("");
			*/
		}
		
		HelperFunctions.outputData(pltFolder + "Apower_Error", "sigma","Error", sigmaNumber, error );

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

}
