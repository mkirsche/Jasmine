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
	static String contigKey = "##contig=<ID=";
	static String filterKey = "##FILTER=<ID=";
	static String fileFormatPrefix = "##fileformat";
	
	ArrayList<String> lines;
	
	HashSet<String> infoFields;
	HashSet<String> formatFields;
	HashSet<String> contigFields;
	HashSet<String> filterFields;
	boolean sawFileFormat;
	
	int lastInfoFieldIndex;
	int lastFormatFieldIndex;
	
	VcfHeader()
	{
		lines = new ArrayList<String>();
		infoFields = new HashSet<String>();
		formatFields = new HashSet<String>();
		contigFields = new HashSet<String>();
		filterFields = new HashSet<String>();
		sawFileFormat = false;
		lastInfoFieldIndex = -1;
		lastFormatFieldIndex = -1;
	}
	
	void print(PrintWriter out)
	{
		for(int i = 0; i<lines.size(); i++)
		{
			out.println(lines.get(i));
		}
	}
	
	void addLine(String line)
	{
		if(line.startsWith(fileFormatPrefix))
		{
			if(sawFileFormat)
			{
				return;
			}
			sawFileFormat = true;
			lines.add(0, line);
			if(lastInfoFieldIndex >= 0)
			{
				lastInfoFieldIndex++;
			}
			if(lastFormatFieldIndex >= 0)
			{
				lastFormatFieldIndex++;
			}
			return;
		}
		
		if(line.startsWith(infoKey))
		{
			String restOfLine = line.substring(infoKey.length());
			int endIdx = restOfLine.indexOf(',');
			String id = restOfLine.substring(0, endIdx);
			
			if(infoFields.contains(id))
			{
				return;
			}
			
			lines.add(line);
			lastInfoFieldIndex = lines.size() - 1;
			lastFormatFieldIndex = lines.size() - 1;
			infoFields.add(id);
			return;
		}
		
		else if(line.startsWith(formatKey))
		{
			String restOfLine = line.substring(formatKey.length());
			int endIdx = restOfLine.indexOf(',');
			String id = restOfLine.substring(0, endIdx);
			
			if(formatFields.contains(id))
			{
				return;
			}
			
			lines.add(line);
			lastFormatFieldIndex = lines.size() - 1;
			formatFields.add(id);
			return;
		}
		
		else if(line.startsWith(contigKey))
		{
			String restOfLine = line.substring(contigKey.length());
			int endIdx = restOfLine.indexOf(',');
			if(endIdx == -1)
			{
				endIdx = restOfLine.indexOf('>');
			}
			String id = restOfLine.substring(0, endIdx);
			
			if(contigFields.contains(id))
			{
				return;
			}
			
			contigFields.add(id);
			lines.add(line);
			return;
		}
		
		else if(line.startsWith(filterKey))
		{
			String restOfLine = line.substring(filterKey.length());
			int endIdx = restOfLine.indexOf(',');
			if(endIdx == -1)
			{
				endIdx = restOfLine.indexOf('>');
			}
			String id = restOfLine.substring(0, endIdx);
			
			if(filterFields.contains(id))
			{
				return;
			}
			
			filterFields.add(id);
			lines.add(line);
			return;
		}
		
		else
		{
			lines.add(line);
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
