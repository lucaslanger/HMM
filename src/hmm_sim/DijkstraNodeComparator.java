package hmm_sim;

import java.util.Comparator;

public class DijkstraNodeComparator implements Comparator<DijkstraNode>{
	public int compare(DijkstraNode d1, DijkstraNode d2) {
		if ( d1.getLengthToNode() < d2.getLengthToNode() ){
			return -1;
		}
		else if(d1.getLengthToNode() == d2.getLengthToNode()){
			return 0;
		}
		else{
			return 1;
		}
		
	}
}
