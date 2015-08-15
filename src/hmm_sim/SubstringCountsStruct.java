package hmm_sim;

import java.util.HashMap;

public class SubstringCountsStruct {
	
	private HashMap<SequenceOfSymbols, HashMap<SequenceOfSymbols, Integer> > counts;

	public SubstringCountsStruct() {
		super();
		this.counts = new HashMap<SequenceOfSymbols, HashMap<SequenceOfSymbols,Integer>>();
	};
	
	public void updateKeyPair(SequenceOfSymbols s1, SequenceOfSymbols s2, int c ){
		if (counts.containsKey(s1)){
			HashMap<SequenceOfSymbols, Integer> h = counts.get(s1);
			h.put(s2, h.get(s2) + c);
			counts.put(s1, h);
		}
		else{
			HashMap<SequenceOfSymbols ,Integer> h = new HashMap<SequenceOfSymbols, Integer>();
			h.put(s2, c);
			counts.put(s1, h);
		}
	}
	
	public SymbolCounts getTotalCounts(){
		SymbolCounts hReturn = new SymbolCounts();
		for (SequenceOfSymbols s : counts.keySet()) {
			for (SequenceOfSymbols seq : counts.get(s).keySet()) {
				hReturn.updateFrequency(seq, counts.get(s).get(seq));
			}
		}
		return hReturn;
		
	}

}
