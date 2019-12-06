/*
 * The header of a VCF file, including a list of INFO fields
 * The main purpose of this is to manage INFO field description lines and avoid duplicates
 */

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashSet;

public class VcfHeader {
	
	static String infoKey = "##INFO=<ID=";
	static String formatKey = "##FORMAT=<ID=";
	
	// List of all lines in the header
	ArrayList<String> lines;
	
	// The names of all INFO fields present in the VCF file
	HashSet<String> infoFields;
	HashSet<String> formatFields;
	
	// The index of the last INFO field line, or the line before where the next INFO field should go
	int lastInfoFieldIndex;
	
	// The index of the last FORMAT field line, or the line before where the next FORMAT field should go
	int lastFormatFieldIndex;
	
	// Constructor - just initializes the data structures
	VcfHeader()
	{
		lines = new ArrayList<String>();
		infoFields = new HashSet<String>();
		formatFields = new HashSet<String>();
		lastInfoFieldIndex = -1;
		lastFormatFieldIndex = -1;
	}
	
	/*
	 * Print all lines of the header
	 */
	void print(PrintWriter out)
	{
		for(int i = 0; i<lines.size(); i++)
		{
			out.println(lines.get(i));
		}
	}
	
	/*
	 * Adds a line to the header, updating the list of INFO fields as needed
	 */
	void addLine(String line)
	{
		lines.add(line);
		
		if(line.startsWith(infoKey))
		{
			lastInfoFieldIndex = lines.size() - 1;
			lastFormatFieldIndex = lines.size() - 1;
			String restOfLine = line.substring(infoKey.length());
			int endIdx = restOfLine.indexOf(',');
			infoFields.add(restOfLine.substring(0, endIdx));
		}
		
		else if(line.startsWith(formatKey))
		{
			lastFormatFieldIndex = lines.size() - 1;
			String restOfLine = line.substring(formatKey.length());
			int endIdx = restOfLine.indexOf(',');
			formatFields.add(restOfLine.substring(0, endIdx));
		}
		
		else
		{
			if(lastFormatFieldIndex == -1)
			{
				lastFormatFieldIndex++;
			}
			if(lastInfoFieldIndex == -1)
			{
				lastInfoFieldIndex++;
			}
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
		lastInfoFieldIndex++;
		lastFormatFieldIndex++;
	}
	
	/*
	 * Remove all format fields from the header
	 */
	void resetFormatFields()
	{
		formatFields = new HashSet<String>();
		ArrayList<String> newLines = new ArrayList<String>();
		int oldIndex = lastFormatFieldIndex;
		lastFormatFieldIndex = -1;
		for(int i = 0; i<lines.size(); i++)
		{
			String line = lines.get(i);
			if(line.startsWith(formatKey))
			{
				if(lastFormatFieldIndex == -1)
				{
					lastFormatFieldIndex = i - 1;
				}
				continue;
			}
			newLines.add(line);
		}
		if(lastFormatFieldIndex == -1)
		{
			lastFormatFieldIndex = oldIndex;
		}
		lines = newLines;
	}
	
	/*
	 * Adds the header line for a new format field, if it doesn't already exist
	 */
	void addFormatField(String id, String number, String type, String desc)
	{
		if(formatFields.contains(id))
		{
			return;
		}
		String line = String.format("##FORMAT=<ID=%s,Number=%s,Type=%s,Description=\"%s\">", id, number, type, desc);
		formatFields.add(id);
		lines.add(lastFormatFieldIndex + 1, line);
		lastFormatFieldIndex++;
	}
}
