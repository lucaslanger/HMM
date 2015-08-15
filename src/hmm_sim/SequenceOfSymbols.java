package hmm_sim;

import java.io.IOException;
import java.io.Serializable;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;

import Jama.Matrix;
import Jama.SingularValueDecomposition;

public class SequenceOfSymbols implements Comparable<SequenceOfSymbols>, Serializable {
	
	private String sequence;
	
	public SequenceOfSymbols(){
	}
	
	public void initialize(String s){
		this.sequence = s;
	}
	
	public SequenceOfSymbols(String s){
		checkSyntaxValidity(s);
		this.sequence = s;
	}
	
	private void checkSyntaxValidity(String s){
		if (s.equals("")){}
		else{
			try{
				String[] streaks = s.split(",");
				for (String string : streaks) {
					String [] sa = string.split(":");
					int symbol = Integer.parseInt(sa[0]);
					int streak = Integer.parseInt(sa[1]);
				}
			}
			catch(Exception e){
				System.out.println("Not a valid string to make a sequence out of");
				e.printStackTrace();
				System.out.println(s);
				System.out.println();
			}
		}
	}
	
	//Should be safe since strings are immutable
	
	public LinkedList<SequenceOfSymbols> getPrefixesFromSequence() { 	//second parameter should be an emptyLinkedlist
		LinkedList<SequenceOfSymbols> L = new LinkedList<SequenceOfSymbols>();
		SequenceOfSymbols s = this;
		while(s.getSequence().equals("") == false){
			L.add( s );
			
			SequenceOfSymbols lastStreak = s.getLastStreak();
			int streak = lastStreak.getStreakFromString();
			String symbol = lastStreak.getSymbolFromString();
			
			if(streak > 1){
				String t = symbol + ":" + Integer.toString(streak-1);
				SequenceOfSymbols seq = new SequenceOfSymbols( s.substring(0, s.rawStringLength() - lastStreak.rawStringLength()) + t );
				s = seq;
			}
			else{	
				s = s.substring(0, s.rawStringLength() - lastStreak.rawStringLength() );		// -1 to Take care of the comma
				if (s.getSequence().equals("") == false && s.getSequence().charAt(s.rawStringLength() - 1) == ','){
					s = s.substring(0, s.rawStringLength() -1 );
				}
			}
			
		}
		L.add( new SequenceOfSymbols("") );

		//System.out.println("Please make sure that sequence did not mutate - From: SequenceOfSymbols.getPrefixesFromSequence");
		//System.out.println(this.sequence);
		
		return L;
		
	}
	
	public static SequenceOfSymbols concatenateSymbols(SequenceOfSymbols s1, SequenceOfSymbols s2){
		if (s1.getSequence().equals( "") ){
			return s2;
		}
		else if (s2.getSequence().equals("")){
			return s1;
		}
		
		SequenceOfSymbols lastStreakOfS1 = s1.getLastStreak();
		SequenceOfSymbols firstStreakOfS2 = s2.getFirstStreak();
		String firstsymbol = lastStreakOfS1.getSymbolFromString();
		String secondsymbol = firstStreakOfS2.getSymbolFromString();
	
		if (firstsymbol.equals(secondsymbol) ){
			int newStreak = firstStreakOfS2.getStreakFromString() + lastStreakOfS1.getStreakFromString();
			String mid = firstsymbol + ":" + Integer.toString(newStreak);
			
			if (s2.rawStringLength() > firstStreakOfS2.rawStringLength()){
				String returnString =  s1.substring(0, s1.rawStringLength() - lastStreakOfS1.rawStringLength()) + mid +  "," + s2.substring(firstStreakOfS2.rawStringLength() + 1, s2.rawStringLength());
				SequenceOfSymbols s = new SequenceOfSymbols(returnString);
				return s;
			}
			
			else{
				String returnString =  s1.substring(0, s1.rawStringLength() - lastStreakOfS1.rawStringLength()) + mid;
				SequenceOfSymbols s = new SequenceOfSymbols(returnString);

				return s;
			}
		}
		else{
			return new SequenceOfSymbols( s1.getSequence() + ","  + s2.getSequence() );
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
		//System.out.println("Weird input or bad behavior");
		//System.out.println(this.sequence + "\n");
		return -1;
	}

	public String getSymbolFromString(){
		
		String streak = "";
		for (int i = 0; i < this.sequence.length(); i++) {
			char c = this.sequence.charAt(i);
			if (c == ':'){
				Integer.parseInt(streak);
				return streak;
			}
			else{
				streak = c + streak;
			}
		}
		if (streak == ""){return streak;}
		
		//System.out.println("Weird input or bad behavior getting symbol");
		//System.out.println(this.sequence);
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
	
	
	public LinkedList<SequenceOfSymbols> getSuffixesOfSequence(){
		SequenceOfSymbols t = this;
		
		LinkedList<SequenceOfSymbols> L = new LinkedList<SequenceOfSymbols>();
		while(t.getSequence().equals("") == false){
			L.add(t);
	
			SequenceOfSymbols firstStreak = t.getFirstStreak();
			int streak = firstStreak.getStreakFromString();
			String symbol = firstStreak.getSymbolFromString();

			String n;
			if (streak > 1){
				
				if (t.rawStringLength() > firstStreak.rawStringLength()){ 
					n = symbol + ":" + Integer.toString(streak-1) + "," + t.substring(firstStreak.rawStringLength()+1, t.rawStringLength() ).getSequence();
					t = new SequenceOfSymbols( n );
				}
				else{
					n = symbol + ":" + Integer.toString(streak-1);
					t = new SequenceOfSymbols( n );
				}
			}
			else{
				if (t.rawStringLength() > firstStreak.rawStringLength() ){	//Delete comma
					n = t.substring(firstStreak.rawStringLength()+1, t.rawStringLength()).getSequence();
					t = new SequenceOfSymbols( n );
				}
				else{
					n = t.substring(firstStreak.rawStringLength(), t.rawStringLength()).getSequence();
					t = new SequenceOfSymbols( n );
					}
				}
		}
		L.add( new SequenceOfSymbols("") );
		
		return L;
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
	
	public int sequenceLength(){
		if (this.sequence.equals("")){
			return 0;
		}
		int c = 0;
		String[] t =  this.sequence.split(",");
		for (String s : t) {
			String l = s.split(":")[1];
			int i = Integer.parseInt(l);
			c += i;
		}
		return c;
	}
	
	public int rawStringLength(){
		return sequence.length();
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
	
	public String getRawSequence(){
		String r = "";
		SequenceOfSymbols t = this;

		while (t.getSequence().equals("") == false){
			SequenceOfSymbols firstStreak = t.getFirstStreak();
			String symbol = firstStreak.getSymbolFromString();
			int streak = firstStreak.getStreakFromString();
			for (int i = 0; i < streak; i++) {
				r = r + symbol;
			}
			if ( t.rawStringLength() > firstStreak.rawStringLength()){
				t = t.substring(firstStreak.rawStringLength()+1, t.rawStringLength() ) ;
			}
			else{
				t = t.substring(firstStreak.rawStringLength(), t.rawStringLength() ); 
			}
		}

		return r;
	}
	
	public static int[] getRawSequencesRowVector(SequenceOfSymbols[] sequences){
		int[] t = new int[sequences.length];
		int i = 0;
		for (SequenceOfSymbols s: sequences){
			String rawSeq = s.getRawSequence();
			if ( rawSeq.equals("") ){
				t[i] = 0;
			}
			else{
				t[i] = Integer.parseInt( rawSeq );
			}
			i++;
		}		
		return t;
	}
	
	@Override
	public boolean equals(Object o){
		//System.out.println("Comparing");
		SequenceOfSymbols s = (SequenceOfSymbols) o;
		if (s.getSequence().equals( this.getSequence() )){
			return true;
		}
		else{
			return false;
		}
	}
	
	public static SequenceOfSymbols fullStringToCompressed(String s){
		if (s.equals("") ){
			return new SequenceOfSymbols("");
		}
		
		int counter = 1;
		String construct = "";
		char prevchar = s.charAt(0);
		for (int i = 1; i < s.length(); i++) {
			if( s.charAt(i) == prevchar ){
				counter++;
			}
			else{
				construct.concat( multiplyChar(counter, prevchar) + "," );
				prevchar = s.charAt(i);
				counter = 1;
			}
			
		}
		construct = construct.concat(multiplyChar(counter, prevchar));
		SequenceOfSymbols rSeq = new SequenceOfSymbols(construct);
		return rSeq;
	}
	
	private static String multiplyChar(int i, char c){
		return c + ":" + i;
	}
	
	private synchronized void writeObject(java.io.ObjectOutputStream stream) throws java.io.IOException{
		stream.writeObject(this.sequence);

	}
	
	private void readObject(java.io.ObjectInputStream in) throws IOException, ClassNotFoundException{
		String sequence = (String) in.readObject();
		this.initialize(sequence);
	}

	public static void printArray(SequenceOfSymbols[] seqs) {
		for (SequenceOfSymbols sequenceOfSymbols : seqs) {
			System.out.print(sequenceOfSymbols.toString() + ", ");
		}
		System.out.println();
	}

	public SymbolCounts getSubstringCounts() {
		SymbolCounts sc = new SymbolCounts();
		String rawSequence = this.getRawSequence();
		for (int i = 0; i < rawSequence.length(); i++) {
			for (int j = i; j < rawSequence.length(); j++) {
				String s = rawSequence.substring(i, j);
				SequenceOfSymbols sCompressed = fullStringToCompressed(s);
				sc.updateFrequency(sCompressed, 1);
			}
		}
		return sc;
	}
	/*
	public SymbolCounts getSubstringCountsNoIntersection() {
	
		SymbolCounts sc = new SymbolCounts();
		for (int i = 0; i < this.getSequence().length()-1; i++) {
			for (int j = i+1; j < this.getSequence().length(); j++) {
				String s = this.getSequence().substring(i, j);
				SequenceOfSymbols seq = new SequenceOfSymbols(s);
				if (sc.getSymbolToFrequency().containsKey(s) == false){
					int lsc = findSubStringCopies(this.getSequence().substring(0, i), s);
					int rsc = findSubStringCopies(this.getSequence().substring(j, this.getSequence().length()), s);		
					sc.updateFrequency( new SequenceOfSymbols(s), lsc + rsc);
				}
			}
		}
	}
	
	public int findSubStringCopies(String s, String searchFor){
		
	}
	*/
	


}
