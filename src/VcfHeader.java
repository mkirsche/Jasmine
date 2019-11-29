/*
 * The header of a VCF file, including a list of INFO fields
 * The main purpose of this is to manage INFO field description lines and avoid duplicates
 */

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

public class VcfHeader {
	
	// List of all lines in the header
	ArrayList<String> lines;
	
	// The names of all INFO fields present in the VCF file
	HashSet<String> infoFields;
	
	int lastInfoFieldIndex;
	
	// Constructor - just initializes the data structures
	VcfHeader()
	{
		lines = new ArrayList<String>();
		infoFields = new HashSet<String>();
		lastInfoFieldIndex = -1;
	}
	
	/*
	 * Print all lines of the header
	 */
	void print(PrintWriter out)
	{
		for(String s : lines)
		{
			out.println(s);
		}
	}
	
	/*
	 * Adds a line to the header, updating the list of INFO fields as needed
	 */
	void addLine(String line)
	{
		lines.add(line);
		
		String infoKey = "##INFO=<ID=";
		if(line.startsWith(infoKey))
		{
			lastInfoFieldIndex = lines.size() - 1;
			String restOfLine = line.substring(infoKey.length());
			int endIdx = restOfLine.indexOf(',');
			infoFields.add(restOfLine.substring(0, endIdx));
		}
	}
	
	/*
	 * Adds an INFO field, adding the header line if it's not already present
	 */
	void addInfoField(String id, String number, String type, String desc)
	{
		if(infoFields.contains(id))
		{
			return;
		}
		String line = String.format("##INFO=<ID=%s,Number=%s,Type=%s,Description=\"%s\">", id, number, type, desc);
		infoFields.add(id);
		lines.add(lastInfoFieldIndex + 1, line);
	}
}
