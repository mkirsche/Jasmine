import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class NormalizeTypes {
	static String inputFile = "";
	static String outputFile = "";
	public static void main(String[] args) throws Exception
	{
		if(args.length != 2)
		{
			System.out.println("Usage: java NormalizeTypes input_vcf output_vcf");
			return;
		}
		else
		{
			inputFile = args[0];
			outputFile = args[1];
			Settings.CHR_NAME_MAP = new ChrNameNormalization();
			convertFile(inputFile, outputFile);
		}		
	}
	
	/*
	 * Convert types in inputFiles to their normalized types and outputs them to a new file
	 */
	static void convertFile(String inputFile, String outputFile) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
				
		PrintWriter out = new PrintWriter(new File(outputFile));
		
		VcfHeader header = new VcfHeader();
		ArrayList<VcfEntry> entries = new ArrayList<VcfEntry>();
				
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith("#"))
			{
				header.addLine(line);
			}
			else if(line.length() > 0)
			{
				VcfEntry ve = VcfEntry.fromLine(line);
				ve.normalizeType();
				entries.add(ve);
			}
		}
				
		header.print(out);
		
		for(VcfEntry ve : entries)
		{
			out.println(ve);
		}
		
		input.close();
		out.close();
	}
}
