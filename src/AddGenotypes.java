/*
 * Adds genotype information to a merged VCF file based on the genotypes of the original variants
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;

public class AddGenotypes {
	
	public static void main(String[] args) throws Exception
	{
		if(args.length == 3)
		{
			String inputFile = args[0];
			String fileList = args[1];
			String outFile = args[2];
			addGenotypes(inputFile, fileList, outFile);
		}
		else
		{
			System.out.println("Usage: java AddGenotypes <inputfile> <vcflist> <outputfile>");
			return;
		}
	}
	
	/*
	 * To add other FORMAT fields, add their details here and add the logic to initialize them in reformatVariantFormat
	 */
	static String[] newFieldNames = {"GT", "IS", "OT", "OS"};
	static String[] newFieldNums = {"1", "1", "1", "1"};
	static String[] newFieldTypes = {"String", "String", "String", "String", "String", "String"};
	static String[] newFieldDescs = new String[] {
			"The genotype of the variant",
			"Whether or not the variant call was marked as specific due to high read support and length",
			"The original type of the variant",
			"The SUPP_VEC field in the previously merged file, if any"
	};
	
	/*
	 * Adds FORMAT fields, including per-sample genotypes, to the variants in a merged VCF file
	 */
	static void addGenotypes(String inputFile, String fileList, String outputFile) throws Exception
	{
		// FORMAT fields of all per-file variant calls
		ArrayList<FileFormatField> inputFormats = new ArrayList<FileFormatField>();
		
		// The names of the samples present across all input files
		ArrayList<String> allSampleNamesList = new ArrayList<String>();
		
		// Get the FORMAT fields of input VCFs by going through the file list
		Scanner input = new Scanner(new FileInputStream(new File(fileList)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			FileFormatField fileFormats = new FileFormatField(line, true);
			for(String sampleName : fileFormats.sampleNames)
			{
				allSampleNamesList.add(inputFormats.size() + "_" + sampleName);
			}
			inputFormats.add(fileFormats);
		}
		input.close();
		
		// Get the number of samples per file to know how much to skip in samples where a variant is absent
		int[] sampleCounts = new int[inputFormats.size()];
		for(int i = 0; i<inputFormats.size(); i++)
		{
			sampleCounts[i] = inputFormats.get(i).sampleNames.length;
		}
				
		// Put all the sample names in an array
		String[] allSampleNames = new String[allSampleNamesList.size()];
		for(int i = 0; i<allSampleNames.length; i++)
		{
			allSampleNames[i] = allSampleNamesList.get(i);
		}
				
		// Now scan through merged VCF and combine FORMAT fields as needed, printing the updated file at the same time
		input = new Scanner(new FileInputStream(new File(inputFile)));
		PrintWriter out = new PrintWriter(new File(outputFile));
		VcfHeader header = new VcfHeader();
		boolean headerPrinted = false;
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
				// If this is the first variant, update and print the header
				if(!headerPrinted)
				{
					headerPrinted = true;
					header.resetFormatFields();
					for(int i = 0; i<newFieldNames.length; i++)
					{
						header.addFormatField(newFieldNames[i], newFieldNums[i], newFieldTypes[i], newFieldDescs[i]);
					}
					
					// Update last line to include all sample names
					String[] lastHeaderLineTokens = header.lines.get(header.lines.size()-1).split("\t");
					if(lastHeaderLineTokens.length >= 9)
					{
						StringBuilder newLastLine = new StringBuilder("");
						for(int i = 0; i<9; i++)
						{
							newLastLine.append(lastHeaderLineTokens[i]);
							if(i < 8)
							{
								newLastLine.append("\t");
							}
						}
						for(String sampleName : allSampleNames)
						{
							newLastLine.append("\t" + sampleName);
						}
						header.lines.set(header.lines.size() - 1, newLastLine.toString());
					}
					
					header.print(out);
				}
				
				// This is the per-variant merging and printing logic
				VcfEntry entry = new VcfEntry(line);
				String suppVec = entry.getInfo("SUPP_VEC");
				if(suppVec.length() == 0)
				{
					// If there is no support vector field, just leave the entry as-is
					out.println(entry);
				}
				else
				{
					// The list of format fields of all variants merged into this one
					ArrayList<VariantFormatField> toMerge = new ArrayList<VariantFormatField>();
					String[] ids = entry.getInfo("IDLIST").split(",");
					for(int i = 0; i<suppVec.length(); i++)
					{
						if(suppVec.charAt(i) == '1')
						{
							String curId = ids[toMerge.size()];
							
							// Get the index of the line where this variant was within its original VCF file
							int variantIndex = inputFormats.get(i).idToVariantIndex.get(curId);
							
							// Add the variant's format fields to the list to merge
							toMerge.add(inputFormats.get(i).variantFormats.get(variantIndex));
						}
					}
					
					// Merge all format fields together and print the resulting VCF entry
					VariantFormatField merged = merge(toMerge, sampleCounts, suppVec);
					for(int i = 0; i<8; i++)
					{
						out.print(entry.tabTokens[i] + "\t");
					}
					out.println(merged);
				}
			}
		}
		input.close();
		out.close();
	}
	
	/*
	 * Merges the format field of multiple variants which share the same FORMAT string
	 * Creates one variant whose set of samples is the concatenation of the inputs' samples
	 */
	static VariantFormatField merge(ArrayList<VariantFormatField> list, int[] sampleCounts, String suppVec)
	{
		int numSamples = 0;
		for(int count : sampleCounts)
		{
			numSamples += count;
		}
		// Initialize empty format field data structure big enough for all of the samples
		VariantFormatField res = new VariantFormatField(numSamples, newFieldNames);
		
		// Update one field at a time
		for(int i = 0; i<newFieldNames.length; i++)
		{
			String fieldName = newFieldNames[i];
			int sampleIndex = 0;
			int listIndex = 0;
			
			// See which samples were in the support vector and fill values accordingly
			for(int j = 0; j<suppVec.length(); j++)
			{
				// Whether or this sample was in the support vector for the variant
				boolean include = suppVec.charAt(j) == '1';
				for(int k = 0; k<sampleCounts[j]; k++)
				{
					if(include)
					{
						// Use the values from the next VariantFormatField
						VariantFormatField cur = list.get(listIndex);
						String val = cur.getValue(k, fieldName);
						if(val.length() > 0)
						{
							res.sampleFieldValues[sampleIndex][res.getFieldIndex(fieldName)] = val;
						}
						else
						{
							res.sampleFieldValues[sampleIndex][res.getFieldIndex(fieldName)] = ".";
						}
						listIndex++;
					}
					else
					{
						// Fill fields with "NA" but use "./." for genotype
						String val = "NA";
						if(fieldName.equals("GT"))
						{
							val = "./.";
						}
						res.sampleFieldValues[sampleIndex][res.getFieldIndex(fieldName)] = val;
					}
					sampleIndex++;
				}
			}
		}
		
		return res;
		
	}
	
	/*
	 * Reformats a variant's format fields to match what we want
	 */
	static VariantFormatField reformatVariantFormat(VariantFormatField oldVariant, VcfEntry entry) throws Exception
	{
		int numSamples = oldVariant.numSamples();
		VariantFormatField res = new VariantFormatField(numSamples, newFieldNames);
		
		for(int i = 0; i<newFieldNames.length; i++)
		{
			String field = newFieldNames[i];
			for(int j = 0; j<numSamples; j++)
			{
				if(field.equals("GT"))
				{
					String oldGt = oldVariant.getValue(j, "GT");
					if(oldGt.length() > 0)
					{
						res.sampleFieldValues[j][i] = oldGt;
					}
					else
					{
						res.sampleFieldValues[j][i] = "./.";
					}
				}
				else if(field.equals("IS"))
				{
					if(entry.hasInfoField("IS_SPECIFIC"))
					{
						res.sampleFieldValues[j][i] = entry.getInfo("IS_SPECIFIC");
					}
					else
					{
						res.sampleFieldValues[j][i] = ".";
					}
				}
				else if(field.equals("OT"))
				{
					if(entry.hasInfoField("OLDTYPE"))
					{
						res.sampleFieldValues[j][i] = entry.getInfo("OLDTYPE");
					}
					else
					{
						String type = entry.getType();
						if(type.length() > 0)
						{
							res.sampleFieldValues[j][i] = type;
						}
						else
						{
							res.sampleFieldValues[j][i] = ".";
						}
					}
				}
				else if(field.equals("OS"))
				{
					if(entry.hasInfoField("SUPP_VEC"))
					{
						res.sampleFieldValues[j][i] = entry.getInfo("SUPP_VEC");
					}
					else
					{
						res.sampleFieldValues[j][i] = ".";
					}
				}
			}
		}
		
		return res;
	}
	
	/*
	 * The values of FORMAT fields for an entire VCF file, including the sample names in the header
	 */
	static class FileFormatField
	{
		// FORMAT field names and value for each individual variant
		ArrayList<VariantFormatField> variantFormats;
		
		// Names of samples which are present in the file
		String[] sampleNames;
		
		// Map from variant ID to index in variantFormats for fast lookup of particular variants
		HashMap<String, Integer> idToVariantIndex;
		
		// The header of the VCF file
		VcfHeader header;
		
		FileFormatField(String fileName, boolean reformat) throws Exception
		{
			variantFormats = new ArrayList<VariantFormatField>();
			idToVariantIndex = new HashMap<String, Integer>();
			header = new VcfHeader();
			Scanner input = new Scanner(new FileInputStream(new File(fileName)));
			boolean extractedSampleNames = false;
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
					// If this is the first variant, we finished the header, so get sample names from the last header line
					if(!extractedSampleNames)
					{
						extractedSampleNames = true;
												
						// Get the list of sample names from the last header line
						String lastLine = header.lines.get(header.lines.size() - 1);
						String[] tabTokens = lastLine.split("\t");
						
						// Check if there are actually sample names in the header
						if(tabTokens.length > 9)
						{
							sampleNames = new String[tabTokens.length - 9];
							for(int i = 0; i<sampleNames.length; i++)
							{
								sampleNames[i] = tabTokens[i + 9];
							}
						}
						else
						{
							sampleNames = new String[0];
						}
					}
					
					// Add this variant's format fields to the list
					VcfEntry entry = new VcfEntry(line);
					idToVariantIndex.put(entry.getId(), variantFormats.size());
					VariantFormatField vff = new VariantFormatField(line);
					if(reformat)
					{
						vff = reformatVariantFormat(vff, entry);
					}
					
					variantFormats.add(vff);
				}
			}
			input.close();
		}
	}
	
	/*
	 * The FORMAT information for a single variant call
	 * It includes the list of fields as well as their values for all samples
	 */
	static class VariantFormatField
	{
		// The names of FORMAT fields in order
		String[] fieldNames;
		
		// The value within each sample of each field, in the same order as in fieldNames 
		String[][] sampleFieldValues;
		
		/*
		 * Initialize the format fields with all "NA" values
		 */
		VariantFormatField(int numSamples, String[] fieldNames)
		{
			this.fieldNames = fieldNames;
			sampleFieldValues = new String[numSamples][fieldNames.length];
			for(int i = 0; i<numSamples; i++)
			{
				Arrays.fill(sampleFieldValues[i], "NA");
			}
		}
		
		/*
		 * Initialize the format fields from a VCF line
		 */
		VariantFormatField(String line) throws Exception
		{
			VcfEntry entry = new VcfEntry(line);
			if(entry.tabTokens.length > 8)
			{
				sampleFieldValues = new String[entry.tabTokens.length - 9][];
				String formatString = entry.tabTokens[8];
				fieldNames = formatString.split(":");
				for(int i = 0; i<sampleFieldValues.length; i++)
				{
					sampleFieldValues[i] = entry.tabTokens[i + 9].split(":");
				}
			}
			else
			{
				sampleFieldValues = new String[0][];
				fieldNames = new String[] {};
			}
		}
		
		/*
		 * Gets the number of samples in the VCF this variant came from
		 */
		int numSamples()
		{
			return sampleFieldValues.length;
		}
		
		/*
		 * Gets the position of a given field in the FORMAT string, or -1 if it's not present
		 */
		int getFieldIndex(String field)
		{
			for(int i = 0; i<fieldNames.length; i++)
			{
				if(fieldNames[i].equals(field))
				{
					return i;
				}
			}
			return -1;
		}
		
		/*
		 * Gets the value of a particular field in the given sample, or "" if it's not present
		 */
		String getValue(int sampleIndex, String field)
		{
			int fieldIndex = getFieldIndex(field);
			if(fieldIndex == -1)
			{
				return "";
			}
			return sampleFieldValues[sampleIndex][fieldIndex];
		}
		
		/*
		 * Gets a VCF-format, tab-separated representation of the FORMAT string plus per-sample genotypes
		 */
		public String toString()
		{
			if(fieldNames.length == 0)
			{
				return "";
			}
			
			StringBuilder res = new StringBuilder("");
			
			// First token is the FORMAT string, with field names separated by ":"
			for(int i = 0; i<fieldNames.length; i++)
			{
				res.append(fieldNames[i]);
				if(i < fieldNames.length - 1)
				{
					res.append(":");
				}
			}
			res.append("\t");
			
			// Field values with samples separated by tabs and values within each sample separated by colons
			for(int i = 0; i<sampleFieldValues.length; i++)
			{
				for(int j = 0; j<sampleFieldValues[i].length; j++)
				{
					res.append(sampleFieldValues[i][j]);
					if(j < sampleFieldValues[i].length - 1)
					{
						res.append(":");
					}
				}
				if(i < sampleFieldValues.length - 1)
				{
					res.append("\t");
				}
			}
			return res.toString();
		}
	}
}
