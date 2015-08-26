package hmm_sim;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;

public class OutputData {

	
	public static void outputData(String filename, String xaxisLabel, String yaxisLabel, double[][] xaxis, double[][] yaxis){
		System.out.println("USING OUTDATA VERSION OF DATA WRITING");
		
		try {
			PrintWriter writer = new PrintWriter(filename, "UTF-8");
			
			writer.println(xaxisLabel + "," + yaxisLabel);
			writer.println("No internal comment");
			writer.println("No title");
			
			StringBuilder line = new StringBuilder();
			for (int j = 0; j < xaxis[0].length; j++) {	
				for (int i = 0; i < xaxis.length; i++) {
					line.append(xaxis[i][j] + ",");
					line.append(yaxis[i][j] + " ");
				}
				writer.println( line );
				line.setLength(0);
			} 
			
			writer.close();
			
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			e.printStackTrace();
		}

	}

	

	public static void outputData(String filename, String xaxisLabel, String yaxisLabel, double[][] xaxis, double[][] yaxis, String plotTitle,  String internalComment){
		try {
			PrintWriter writer = new PrintWriter(filename, "UTF-8");
			
			writer.println(xaxisLabel + "," + yaxisLabel);
			writer.println(internalComment);
			writer.println(plotTitle);
			
			StringBuilder line = new StringBuilder();
			for (int j = 0; j < xaxis[0].length; j++) {	
				for (int i = 0; i < xaxis.length; i++) {
					line.append(xaxis[i][j] + ",");
					line.append(yaxis[i][j] + " ");
				}
				writer.println( line );
				line.setLength(0);
			} 
			
			writer.close();
			
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			System.out.println("Something messed up when outputting data");
			e.printStackTrace();
		}

	}
	
	
	public static void outputData(String filename, String xaxisLabel, String yaxisLabel, double[][] xaxis, double[][] yaxis, double[][] spreads, String plotTitle,  String internalComment){
		try {
			PrintWriter writer = new PrintWriter(filename, "UTF-8");
			
			writer.println(xaxisLabel + "," + yaxisLabel);
			writer.println(internalComment);
			writer.println(plotTitle);
			
			StringBuilder line = new StringBuilder();
			for (int j = 0; j < xaxis[0].length; j++) {	
				for (int i = 0; i < xaxis.length; i++) {
					line.append(xaxis[i][j] + ",");
					line.append(yaxis[i][j] + ",");
					line.append(spreads[i][j] + " ");
				}
				writer.println( line );
				line.setLength(0);
			} 
			
			writer.close();
			
		} catch (FileNotFoundException | UnsupportedEncodingException e) {
			System.out.println("Something messed up when outputting data");
			e.printStackTrace();
		}

	}
}
