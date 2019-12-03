/*
 * Adds genotype information to a merged VCF file based on the genotypes of the original variants
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
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
	 * Adds per-sample genotypes to all of the variants in a merged sample
	 */
	@SuppressWarnings("unchecked")
	static void addGenotypes(String inputFile, String fileList, String outputFile) throws Exception
	{
		// The number of VCF files used in merging
		int n = VariantInput.countFiles(fileList);
		
		// Get every field's name and values from the input VCFs
		ArrayList<GTField>[] fields = new ArrayList[n];
		for(int i = 0; i<n; i++) fields[i] = new ArrayList<GTField>();
		Scanner fileListInput = new Scanner(new FileInputStream(new File(fileList)));
		int fileIndex = 0;
		while(fileListInput.hasNext())
		{
			String vcfFile = fileListInput.nextLine();
			if(vcfFile.length() == 0)
			{
				continue;
			}
			
			fields[fileIndex] = getSampleGenotypeFields(vcfFile);
			fileIndex++;
		}
		fileListInput.close();
		
		// All of the names of GT fields as a single list for use in the last line of the header
		String[] sampleNames = getAllGtFieldNames(fields);
		
		// Scan through the merged VCF and print modified lines to the output VCF
		Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
		PrintWriter out = new PrintWriter(new File(outputFile));
		VcfHeader header = new VcfHeader();
		boolean printedHeader = false;
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
				if(!printedHeader)
				{
					header.print(out, sampleNames);
					printedHeader = true;
				}
				VcfEntry entry = new VcfEntry(line);
				String suppVec = entry.getInfo("SUPP_VEC");
				String[] ids = entry.getInfo("IDLIST").split(",");
				int idIdx = 0;
				for(int i = 0; i<9; i++)
				{
					out.print(entry.tabTokens[i]);
					if(i < 8)
					{
						out.print("\t");
					}
				}
				
				// Go through each sample and print the genotype of this variant in that sample
				for(int i = 0; i<suppVec.length(); i++)
				{
					boolean inSample = suppVec.charAt(i) == '1';
					for(int j = 0; j<fields[i].size(); j++)
					{
						if(inSample)
						{
							// Print the genotype values of this field from the corresponding sample 
							out.print("\t" + fields[i].get(j).idToVal.get(ids[idIdx]));
						}
						else
						{
							// Not in this sample, so print NA for every value
							String[] formatNames = entry.tabTokens[8].split(":");
							for(int c = 0; c<formatNames.length; c++)
							{
								if(c == 0) out.print("\t");
								else out.print(":");
								if(formatNames[c].equalsIgnoreCase("GT"))
								{
									out.print("./.");
								}
								else
								{
									out.print("NA");
								}
							}
						}
					}
					if(inSample)
					{
						idIdx++;
					}
				}
				out.println();
			}
		}
		
		input.close();
		out.close();
	}
	
	/*
	 * Get the genotype fields from a VCF file, along with each field's value for every variant
	 */
	static ArrayList<GTField> getSampleGenotypeFields(String vcfFile) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(vcfFile)));
		
		// Get the name for GT columns from this file based on the filename
		String curColName = StringUtils.fileBaseName(vcfFile);
		if(curColName.endsWith(".vcf") || curColName.endsWith(".VCF"))
		{
			curColName = curColName.substring(0, curColName.length() - 4);
		}
		
		// Some bookkeeping variables for getting and parsing the last header line
		String lastLine = null;
		boolean processedHeader = false;
		
		// The number of GT fields in the file - will be initialized once the full header is parsed
		int numFields = 0;
		
		// List of all GT fields and their values for every variant
		ArrayList<GTField> curFields = new ArrayList<GTField>();
		
		// Scan through the file
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			
			// Keep track of the most recent header line so we can work with the last header line once we know the header is done
			if(line.startsWith("#"))
			{
				lastLine = line;
			}
			else
			{
				
				// If this is the first non-header line, use the last header line to get the names of all GT fields
				if(!processedHeader)
				{
					String[] lastTokens = lastLine.split("\t");
					numFields = lastTokens.length - 9;
					if(numFields <= 0)
					{
						break;
					}
					else
					{
						processedHeader = true;
						if(numFields == 1)
						{
							// Name the only GT field based on the file's basename
							curFields.add(new GTField(curColName));
						}
						else
						{
							// Use the original field names to differentiate between them 
							for(int i = 0; i<numFields; i++)
							{
								curFields.add(new GTField(curColName +  "_" + lastTokens[i+9]));
							}
						}
					}
					processedHeader = true;
				}
				
				// Process this variant and add its value for each GT field to the map
				VcfEntry entry = new VcfEntry(line);
				for(int i = 0; i<numFields; i++)
				{
					String id = entry.getId();
					String val = entry.tabTokens[i+9];
					curFields.get(i).addVariant(id, val);
				}
			}
		}
		input.close();
		
		return curFields;
	}
	
	/*
	 * Gets the name of all GT fields across samples as a single array
	 */
	static String[] getAllGtFieldNames(ArrayList<GTField>[] all)
	{
		ArrayList<String> resList = new ArrayList<String>();
		for(ArrayList<GTField> sampleFields : all)
		{
			for(GTField field : sampleFields)
			{
				resList.add(field.name);
			}
		}
		String[] res = new String[resList.size()];
		for(int i = 0; i<res.length; i++)
		{
			res[i] = resList.get(i);
		}
		return res;
	}
	
	/*
	 * The values within a single VCF file of one GT field
	 */
	static class GTField
	{
		String name;
		HashMap<String, String> idToVal;
		GTField(String name)
		{
			this.name = name;
			idToVal = new HashMap<String, String>();
		}
		void addVariant(String id, String gtVal)
		{
			idToVal.put(id, gtVal);
		}
	}
}
