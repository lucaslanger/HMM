package hmm_sim;

import java.util.Comparator;

public class SymbolProbabilityPairComparator implements Comparator<SymbolProbabilityPair>{

	
	public int compare(SymbolProbabilityPair p1, SymbolProbabilityPair p2) {
		if ( p1.getProbability() < p2.getProbability() ){
			return -1;
		}
		else if( p1.getProbability() == p2.getProbability()){
			return 0;
		}
		else{
			return 1;
		}
		
	}
}
