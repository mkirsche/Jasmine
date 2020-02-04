/*
 * Program for marking calls which fall in a specific callset, given the sensitive callset
 * Takes parameters for read support and SV length required for a variant to be specific
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Scanner;

public class MarkSpecificCalls {
public static void main(String[] args) throws Exception
{
	String fn = "";
	String ofn = "";
	int minReadSupport = 0;
	int minLength = 0;
	if(args.length == 4)
	{
		fn = args[0];
		ofn = args[1];
		minReadSupport = Integer.parseInt(args[2]);
		minLength = Integer.parseInt(args[3]);
		convertFile(fn, ofn, minReadSupport, minLength);
	}
	else
	{
		System.out.println("Usage: java MarkSpecificCalls vcffile outfile minreadsupport minlength");
		return;
	}	
}

/*
 * Marks specific calls in inputFile and outputs updated VCF to outputFile
 * A variant will be specific if its number of supporting reads is at least
 * minReadSupport (or unspecified) and its length (absolute value) is at least minLength
 */
static void convertFile(String inputFile, String outputFile, int minReadSupport, int minLength) throws Exception
{
	Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
	PrintWriter out = new PrintWriter(new File(outputFile));
	
	VcfHeader header = new VcfHeader();
	ArrayList<VcfEntry> entries = new ArrayList<VcfEntry>();
	
	while(input.hasNext())
	{
		String line = input.nextLine();
		if(line.length() == 0)
		{
			continue;
		}
		if(line.startsWith("#"))
		{
			header.addLine(line);
		}
		else
		{
			VcfEntry entry = VcfEntry.fromLine(line);
			boolean inSpecific = false;
			int readSupport = entry.getReadSupport();
			
			boolean longEnough = entry.getType().equals("TRA") || entry.getType().equals("BND") || Math.abs(entry.getLength()) >= minLength || entry.getLength() == 0;
			
			if(readSupport >= minReadSupport && longEnough)
			{
				inSpecific = true;
			}
			
			entry.setInfo("IN_SPECIFIC", inSpecific ? "1" : "0");
			entries.add(entry);
		}
	}
	
	header.addInfoField("IN_SPECIFIC", "1", "String", "Whether or not a variant has enough read support and length to be specific");
	header.print(out);
	
	for(VcfEntry entry : entries)
	{
		out.println(entry);
	}
	
	input.close();
	out.close();
}
}
