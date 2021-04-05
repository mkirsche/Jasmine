/*
 * Converts insertions which were originally duplications back to their original SV calls
 * Usage: java InsertionsToDuplications input_vcf output_vcf
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class InsertionsToDuplications {
	static String inputFile = "";
	static String outputFile = "";
	public static void main(String[] args) throws Exception
	{
		if(args.length != 2)
		{
			System.out.println("Usage: java InsertionsToDuplications input_vcf output_vcf");
			return;
		}
		else
		{
			inputFile = args[0];
			outputFile = args[1];
			convertFile(inputFile, outputFile);
		}		
	}
	
	/*
	 * Convert any insertions which have OLDTYPE marked as DUP back to duplications
	 */
	static void convertFile(String inputFile, String outputFile) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
		
		PrintWriter out = new PrintWriter(new File(outputFile));
		
		VcfHeader header = new VcfHeader();
		
		int countDup = 0;
		
		boolean headerPrinted = false;
		int totalEntries = 0;
		
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith("#"))
			{
				header.addLine(line);
				continue;
			}
			
			if(!headerPrinted)
			{
				header.addInfoField("REFINEDALT", "1", "String", "For duplications which were changed to insertions and refined, the refined ALT sequence");
				header.addInfoField("STRANDS", "1", "String", "");
				header.print(out);
				headerPrinted = true;
			}
			
			VcfEntry ve = new VcfEntry(line);
			
			totalEntries++;
			
			if(line.contains("OLDTYPE=DUP") && ve.getType().equals("INS"))
			{
				countDup++;
					
				long start = ve.getPos();
				int length = ve.getLength();
				long nstart = start - length + 1, nend = nstart + length;
				String refinedAlt = ve.getAlt();
				ve.setPos(nstart);
				ve.setInfo("END", nend+"");
				ve.setType("DUP");
				ve.setInfo("REFINEDALT", refinedAlt);
				ve.setInfo("STRANDS", "-+");
				ve.setRef(".");
				ve.setAlt("<DUP>");
				out.println(ve);
			}
			else
			{
				ve.setInfo("REFINEDALT", ".");
				out.println(ve);
			}
		}
		
		System.out.println("Number of insertions converted back to duplications: " + countDup + " out of " + totalEntries + " total variants");
		
		input.close();
		out.close();
	}
}
