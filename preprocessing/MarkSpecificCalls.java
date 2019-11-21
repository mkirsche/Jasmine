/*
 * Program for marking calls which fall in a specific callset, given the sensitive callset
 * Takes parameters for read support and SV length required for the specific callset
 */
import java.util.*;
import java.io.*;
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
	}
	else
	{
		System.out.println("Usage: java MarkSpecificCalls vcffile outfile minreadsupport minlength");
		return;
	}
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	PrintWriter out = new PrintWriter(new File(ofn));
	
	ArrayList<String> headerLines = new ArrayList<String>();
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
			headerLines.add(line);
		}
		else
		{
			VcfEntry entry = new VcfEntry(line);
			boolean inSpecific = false;
			int readSupport = 0;
			String[] rnamesField = entry.getRnames();
			if(rnamesField.length > 0)
			{
				readSupport = rnamesField.length;
			}
			
			if(readSupport >= minReadSupport && Math.abs(entry.getLength()) >= minLength) 
			{
				inSpecific = true;
			}
			
			entry.setInfo("IN_SPECIFIC", inSpecific ? "1" : "0");
			entries.add(entry);
		}
	}
	
	for(String s : headerLines)
	{
		out.println(s);
	}
	
	out.println("##INFO=<ID=IN_SPECIFIC,Number=1,Type=String,Description=\"\">");
	
	for(VcfEntry entry : entries)
	{
		out.println(entry);
	}
	
	input.close();
	out.close();
}
}
