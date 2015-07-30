package hmm_sim;

import java.util.LinkedList;

public class SequenceOfSymbols implements Comparable<SequenceOfSymbols> {
	
	private String sequence;
	
	public SequenceOfSymbols(String s){
		this.sequence = s;
	}
	
	//Should be safe since strings are immutable
	
	public LinkedList<SequenceOfSymbols> getPrefixesFromSequence() { 	//second parameter should be an emptyLinkedlist
		LinkedList<SequenceOfSymbols> currentList = new LinkedList<SequenceOfSymbols>();
		SequenceOfSymbols s = this;
		while(s.getSequence().equals("") == false){
			
			currentList.add( s );
			
			SequenceOfSymbols lastStreak = s.getLastStreak();
			int streak = lastStreak.getStreakFromString();
			String symbol = lastStreak.getSymbolFromString();
			
			if(streak > 1){
				String t = symbol + ":" + Integer.toString(streak-1);
				SequenceOfSymbols seq = new SequenceOfSymbols( s.substring(0, s.length() - lastStreak.length()).getSequence() + t );
				s = seq;
			}
			else{	
				s = s.substring(0, s.length() - lastStreak.length() );		// -1 to Take care of the comma
				if (s.getSequence().charAt(s.length()-1) == ','){
					s = s.substring(0, s.length() -1 );
				}
			}
			
		}
		System.out.println("Please make sure that sequence did not mutate");
		System.out.println(this.sequence);
		
		return currentList;
		
	}
	
	public SequenceOfSymbols concatenateSymbols(SequenceOfSymbols s1, SequenceOfSymbols s2){
		SequenceOfSymbols lastStreakOfS1 = s1.getLastStreak();
		SequenceOfSymbols firstStreakOfS2 = s2.getFirstStreak();
		String firstsymbol = lastStreakOfS1.getSymbolFromString();
		String secondsymbol = firstStreakOfS2.getSymbolFromString();
		if (firstsymbol.equals(secondsymbol) ){
			int newStreak = firstStreakOfS2.getStreakFromString() + lastStreakOfS1.getStreakFromString();
			
			String mid = firstsymbol + ":" + Integer.toString(newStreak);
			
			String returnString =  s1.substring(0, s1.length() - lastStreakOfS1.length()) + mid + s2.substring(firstStreakOfS2.length(), s2.length());
			SequenceOfSymbols s = new SequenceOfSymbols(returnString);
			return s;
		}
		else{
			if(s2.getSequence() != ""){
				return new SequenceOfSymbols(s1.getSequence() + "," + s2.getSequence());
			}
			else{
				return s1;
			}
		}
	}
	
	public int getStreakFromString() {
		String streak = "";
		for (int i = this.sequence.length()-1; i >= 0  ; i--) {
			char c = this.sequence.charAt(i);
			if (c == ':'){
				return Integer.parseInt(streak);
			}
			else{
				streak = c + streak;
			}
		}
		System.out.println("Weird input or bad behavior");
		System.out.println(this.sequence + "\n");
		return -1;
	}

	public String getSymbolFromString(){
		
		String streak = "";
		for (int i = 0; i < this.sequence.length(); i++) {
			char c = this.sequence.charAt(i);
			if (c == ':'){
				return streak;
			}
			else{
				streak = c + streak;
			}
		}
		System.out.println("Weird input or bad behavior");
		return null;
	}

	public SequenceOfSymbols getFirstStreak() {
		String r = "";
		for (int i = 0; i < this.sequence.length(); i++) {
			char c = this.sequence.charAt(i);
			if(c == ','){
				return new SequenceOfSymbols(r);
			}
			else{
				r = r + c;
			}
		}
		return new SequenceOfSymbols(r);
	}

	private SequenceOfSymbols getLastStreak() {
		String r = "";
		for (int i = this.sequence.length()-1; i >= 0; i--) {
			char c = this.sequence.charAt(i);
			if (c == ',') {
				return new SequenceOfSymbols(r);
			} 
			else {
				r = c + r;
			}
		}
		return new SequenceOfSymbols(r);
	}
	
	public SequenceOfSymbols substring(int start, int finish){
		return new SequenceOfSymbols(this.sequence.substring(start, finish));
	}
	
	public String getSequence(){
		return this.sequence;
	}
	
	public int length(){
		return sequence.length();
	}
	
	public static SequenceOfSymbols concatenate(SequenceOfSymbols s1, SequenceOfSymbols s2){
		return s2;
	}
	
	public int hashCode(){
		return sequence.hashCode();
	}

	@Override
	public int compareTo(SequenceOfSymbols s2) {
		return this.sequence.compareTo(s2.getSequence());
	}
	
	public String toString(){
		return sequence;
	}


}
