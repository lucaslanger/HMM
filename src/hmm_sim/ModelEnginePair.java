package hmm_sim;

import java.util.HashMap;

public class ModelEnginePair {

	private HashMap<Integer, QueryEngine[][]> anySizeQueryEngines;
	private HashMap<Integer, HankelSVDModel[]> empericalModels;
	
	public ModelEnginePair(HashMap<Integer, QueryEngine[][]> anySizeQueryEngines, HashMap<Integer, HankelSVDModel[]> empericalModels){
		this.anySizeQueryEngines = anySizeQueryEngines;
		this.empericalModels = empericalModels;
	}
	
	public HashMap<Integer, QueryEngine[][]> getAnySizeQueryEngines() {
		return anySizeQueryEngines;
	}

	public HashMap<Integer, HankelSVDModel[]> getEmpericalModels() {
		return empericalModels;
	}
}

