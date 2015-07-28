package hmm_sim;

import java.util.Comparator;

public class SymbolCountPairComparator implements Comparator<SymbolCountPair>{

	
	public int compare(SymbolCountPair p1, SymbolCountPair p2) {
		if ( p1.getCount() < p2.getCount() ){
			return -1;
		}
		else if( p1.getCount() == p2.getCount()){
			return 0;
		}
		else{
			return 1;
		}
		
	}
}
