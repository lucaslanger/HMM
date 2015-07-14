package hmm_sim;

public class DijkstraNode{

	private int id;
	private int lengthToNode;
	
	public DijkstraNode(int id, int lengthToNode) {
		this.id = id;
		this.lengthToNode = lengthToNode;
	}
	
	public int getId() {
		return id;
	}

	public int getLengthToNode() {
		return lengthToNode;
	}
	
	public void setLengthToNode(int l){
		this.lengthToNode = l;
	}

}
