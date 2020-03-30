/*
 * Converts duplications to insertions
 * Usage: java DuplicationsToInsertions input_vcf reference_genome output_vcf
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class DuplicationsToInsertions {
	static String inputFile = "";
	static String genomeFile = "";
	static String outputFile = "";
	public static void main(String[] args) throws Exception
	{
		if(args.length != 3)
		{
			System.out.println("Usage: java DuplicationsToInsertions input_vcf reference_genome output_vcf");
			return;
		}
		else
		{
			inputFile = args[0];
			genomeFile = args[1];
			outputFile = args[2];
			convertFile(inputFile, genomeFile, outputFile);
		}		
	}
	
	/*
	 * Convert duplications in inputFile to insertions and write updated VCF to outputFile.
	 * A genome file is needed to get the insertion sequences based on the duplication position
	 */
	static void convertFile(String inputFile, String genomeFile, String outputFile) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
		
		GenomeQuery gq = new GenomeQuery(genomeFile);
		
		PrintWriter out = new PrintWriter(new File(outputFile));
		
		VcfHeader header = new VcfHeader();
		ArrayList<VcfEntry> entries = new ArrayList<VcfEntry>();
		
		int countDup = 0;
		
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith("#"))
			{
				header.addLine(line);
			}
			else if(line.length() > 0)
			{
				VcfEntry ve = new VcfEntry(line);
				if(ve.getType().equals("DUP") && ve.getLength() < Settings.MAX_DUP_LEN)
				{
					countDup++;
					
					long start = ve.getPos(), end = Long.parseLong(ve.getInfo("END"));
					int length = ve.getLength();
					long nstart = start + length - 1, nend = nstart;

					if(ve.getAlt().equals("<DUP>"))
					{
						String seq = gq.genomeSubstring(ve.getChromosome(), start, end-1);
						
						if(length < 100000)
						{
							ve.setRef(seq.charAt(seq.length()-1)+"");
							ve.setAlt(seq.charAt(seq.length()-1)+seq);
						}
						else
						{
							ve.setRef(".");
							ve.setAlt("<INS>");
						}
						ve.setInfo("END", nend+"");
						ve.setInfo("STRANDS", "+-");
						ve.setPos(nstart);
						ve.setInfo("OLDTYPE", "DUP");
					}
					else
					{
						ve.setInfo("OLDTYPE", "DUP");
					}
					ve.setType("INS");
					
				}
				else
				{
					ve.setInfo("OLDTYPE", ve.getType());
				}
				entries.add(ve);
			}
		}
		
		System.out.println("Number of duplications converted to insertions: " + countDup + " out of " + entries.size() + " total variants");
		
		header.addInfoField("OLDTYPE", "1", "String", "");
		header.addInfoField("STRANDS", "1", "String", "");
		header.print(out);
		
		for(VcfEntry ve : entries)
		{
			out.println(ve);
		}
		
		input.close();
		out.close();
	}
}
