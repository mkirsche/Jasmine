/*
 * Converts insertions which were originally duplications back to their original SV
 * Usage: java InsertionsToDuplications input_vcf output_vcf
 */
import java.util.*;
import java.io.*;
public class InsertionsToDuplications {
	static String inputFile = "/home/mkirsche/entex/enc003_nodups.refined.vcf";
	static String outputFile = "/home/mkirsche/entex/enc003_dups_refined.vcf";
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
		}
		
		convertFile(inputFile, outputFile);
	}
	static void convertFile(String inputFile, String outputFile) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
		
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
			else if(line.contains("OLDTYPE=DUP"))
			{
				VcfEntry ve = new VcfEntry(line);
				countDup++;
					
				long start = ve.getPos();
				int length = ve.getLength();
				long nstart = start - length, nend = nstart;
				String refinedAlt = ve.getAlt();
				ve.setPos(nstart);
				ve.setInfo("END", nend+"");
				ve.setType("DUP");
				ve.setInfo("REFINEDALT", refinedAlt);
				ve.setRef(".");
				ve.setAlt("<DUP>");
				entries.add(ve);
			}
			else
			{
				VcfEntry ve = new VcfEntry(line);
				ve.setInfo("REFINEDALT", ".");
				entries.add(ve);
			}
		}
		
		System.out.println("Number of insertions converted back to duplications: " + countDup + " out of " + entries.size() + " total variants");
		
		header.addInfoField("REFINEDALT", "1", "String", "");
		header.print(out);
						
		for(VcfEntry ve : entries)
		{
			out.println(ve);
		}
		
		input.close();
		out.close();
	}
}
