/*
 * A data structure for holding merged variants to output
 * When merging is performed, the resulting variants use information from multiple files,
 * so some bookkeeping is required to scan through the files one at a time and update all merged variants at once
 */

import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeMap;

public class VariantOutput {
	
	TreeMap<String, VariantGraph> groups;
	
	VariantOutput()
	{
		groups = new TreeMap<String, VariantGraph>();
	}
	
	/*
	 * Adds a graph to the output so that it can be accessed by graph ID easily
	 */
	void addGraph(String graphID, ArrayList<Variant>[] graph, int sampleCount)
	{
		groups.put(graphID,  new VariantGraph(graph, sampleCount));
	}
	
	/*
	 * Given a list of VCF files and merging results, output an updated VCF file
	 */
	public void writeMergedVariants(String fileList, String outFile) throws Exception
	{
		PrintWriter out = new PrintWriter(new File(outFile));
		int sample = 0;
		
		VcfHeader header = new VcfHeader();
		
		boolean printedHeader = false;
		
		ArrayList<String> filenames = PipelineManager.getFilesFromList(fileList);
		for(String filename : filenames)
		{
			Scanner input = new Scanner(new FileInputStream(new File(filename)));
			
			// Iterate over the variants in that file
			while(input.hasNext())
			{
				String line = input.nextLine();
				
				// Ignore empty lines
				if(line.length() == 0)
				{
					continue;
				}
				
				// Print header lines from the first file listed
				else if(line.startsWith("#"))
				{
					if(sample == 0)
					{
						header.addLine(line);
					}
					else
					{
						continue;
					}
				}
				else
				{
					
					// Update the consensus variant in the appropriate graph
					if(sample == 0 && !printedHeader)
					{
						printedHeader = true;
						header.addInfoField("SUPP_VEC", "1", "String", "Vector of supporting samples");
						header.addInfoField("SUPP_VEC_EXT", "1", "String", "Vector of supporting samples, potentially extended across multiple merges");
						header.addInfoField("SUPP", "1", "String", "Number of samples supporting the variant");
						header.addInfoField("SUPP_EXT", "1", "String", "Number of samples supporting the variant, potentially extended across multiple merges");
						header.addInfoField("IDLIST", ".", "String", "Variant IDs of variants merged to make this call");
						header.addInfoField("IDLIST_EXT", ".", "String", "Variant IDs of variants merged, potentially extended across multiple merges");
						header.addInfoField("SVMETHOD", "1", "String", "");
						header.addInfoField("STARTVARIANCE", "1", "String", "Variance of start position for variants merged into this one");
						header.addInfoField("ENDVARIANCE", "1", "String", "Variance of end position for variants merged into this one");
						header.addInfoField("END", "1", "String", "The end position of the variant");
						header.addInfoField("SVLEN", "1", "String", "The length (in bp) of the variant");
						header.print(out);
					}
					VcfEntry entry = new VcfEntry(line);
					String graphID = entry.getGraphID();
					groups.get(graphID).processVariant(entry, sample, out);
				}
			}
			input.close();
			sample++;
		}
		
		out.close();
	}
	
	/*
	 * A graph of variants which are all on the same chromosome, and have the same type and/or strand as specified by the user
	 * This class has login for storing and updating the connected components of the graph as consensus variants
	 */
	class VariantGraph
	{
		// Number of variants in each group
		int[] sizes;
		
		// How many variants in each group have been seen so far
		int[] used;
		
		// The current consensus variant for each group
		VcfEntry[] consensus;
		
		// For each variant ID, the group number it is in
		HashMap<String, Integer> varToGroup;
		
		// For each group, the start, end, and length of its first variant
		HashMap<Integer, Long> firstStart;
		HashMap<Integer, Long> firstEnd;
		HashMap<Integer, Integer> firstLength;
		
		// For each group, the support vector of samples it's in
		String[] supportVectors;
		
		// For each group, the number of sample it's in
		int[] supportCounts;
		
		// The list of variant IDs in each merged variant
		StringBuilder[] idLists;
		
		VariantGraph(ArrayList<Variant>[] groups, int sampleCount)
		{
			int n = groups.length;
			sizes = new int[n];
			used = new int[n];
			consensus = new VcfEntry[n];
			supportVectors = new String[n];
			supportCounts = new int[n];
			idLists = new StringBuilder[n];
			varToGroup = new HashMap<String, Integer>();
			
			// Initialize information about first variant in each group if using it
			if(Settings.KEEP_FIRST_POS)
			{
				firstStart = new HashMap<Integer, Long>();
				firstEnd = new HashMap<Integer, Long>();
				firstLength = new HashMap<Integer, Integer>();
			}
			
			// Scan through groups and map variant IDs to group numbers
			for(int i = 0; i<n; i++)
			{
				sizes[i] = groups[i].size();
				consensus[i] = null;
				idLists[i] = new StringBuilder("");
				char[] suppVec = new char[sampleCount];
				Arrays.fill(suppVec, '0');
				for(int j = 0; j<sizes[i]; j++)
				{
					int sampleID = groups[i].get(j).sample;
					if(suppVec[sampleID] == '0')
					{
						suppVec[sampleID] = '1';
						supportCounts[i]++;
					}
					String idString = groups[i].get(j).id;
					varToGroup.put(idString, i);
				}
				supportVectors[i] = new String(suppVec);
			}
		}
		
		/*
		 * From a VCF line, update the appropriate consensus entry
		 */
		void processVariant(VcfEntry entry, int sample, PrintWriter out) throws Exception
		{
			// This should never happen, but if the variant ID is not in the graph ignore it
			String fullId = VariantInput.fromVcfEntry(entry, sample).id;
			if(!varToGroup.containsKey(fullId))
			{
				return;
			}
			
			int groupNumber = varToGroup.get(fullId);
			
			// Don't even store the components with too little support to be output
			if(supportCounts[groupNumber] < Settings.MIN_SUPPORT)
			{
				return;
			}
			
			// If this is the first variant in the group, initialize the consensus entry
			if(used[groupNumber] == 0)
			{
				consensus[groupNumber] = entry;
				consensus[groupNumber].setId(fullId);
				
				String varId = entry.oldId;
				
				idLists[groupNumber].append(varId);
				consensus[groupNumber].setInfo("END", entry.getEnd() + "");
				consensus[groupNumber].setInfo("SVLEN", entry.getLength() + "");
				consensus[groupNumber].setInfo("STARTVARIANCE", (entry.getPos() * entry.getPos()) + "");
				consensus[groupNumber].setInfo("ENDVARIANCE", (entry.getEnd() * entry.getEnd()) + "");
				
				if(Settings.KEEP_FIRST_POS)
				{
					firstStart.put(groupNumber, entry.getPos());
					firstEnd.put(groupNumber, entry.getEnd());
					firstLength.put(groupNumber, entry.getLength());
				}
				
				if(Settings.INPUTS_MERGED)
				{
					// Set the extended support vector
					// First look for extended suppVec fields, then regular suppVec fields, then just use a single "1"
					String suppVecExt = entry.getInfo("SUPP_VEC_EXT");
					if(suppVecExt.length() == 0)
					{
						suppVecExt = entry.getInfo("SUPP_VEC");
					}
					if(suppVecExt.length() == 0)
					{
						suppVecExt = "1";
					}
					
					// Add zeroes for any absent sample before this
					for(int i = 0; i<sample; i++)
					{
						int previouslyMergedCount = VariantInput.previouslyMergedSamples.get(i);
						for(int j = 0; j<previouslyMergedCount; j++)
						{
							suppVecExt = "0" + suppVecExt;
						}
					}
					consensus[groupNumber].setInfo("SUPP_VEC_EXT", suppVecExt);
					
					// Set the extended ID List
					// First look for extended IDList fields, then regular IDList fields
					String idListExt = entry.getInfo("IDLIST_EXT");
					if(idListExt.length() == 0)
					{
						idListExt = entry.getInfo("IDLIST");
					}
					if(idListExt.length() == 0)
					{
						idListExt = entry.oldId;
					}
	
					consensus[groupNumber].setInfo("IDLIST_EXT", idListExt);
					
					// Update extended support count
					int extendedSupport = 1;
					if(entry.hasInfoField("SUPP_EXT"))
					{
						extendedSupport = Integer.parseInt(entry.getInfo("SUPP_EXT"));
					}
					else if(entry.hasInfoField("SUPP"))
					{
						extendedSupport = Integer.parseInt(entry.getInfo("SUPP"));
					}
					consensus[groupNumber].setInfo("SUPP_EXT", extendedSupport + "");
				}
			}
			
			// Otherwise, update the consensus to include info from this variant
			// For average, don't divide yet to avoid loss of precision
			else
			{
				if(entry.getInfo("OLDTYPE").equals("DUP"))
				{
					consensus[groupNumber].setInfo("OLDTYPE", "DUP");
				}
				
				// Update start
				consensus[groupNumber].setPos(consensus[groupNumber].getPos() + entry.getPos());
				
				// Update end
				long oldEnd = consensus[groupNumber].getEnd();
				long newEnd = entry.getEnd();
				consensus[groupNumber].setInfo("END", (oldEnd + newEnd) + "");
				
				// Update start/end variance (stored as sum of squares of start/end - variance is computed at the end)
				long oldStartVar = Long.parseLong(consensus[groupNumber].getInfo("STARTVARIANCE"));
				long oldEndVar = Long.parseLong(consensus[groupNumber].getInfo("ENDVARIANCE"));
				consensus[groupNumber].setInfo("STARTVARIANCE", (oldStartVar + entry.getPos() * entry.getPos()) + "");
				consensus[groupNumber].setInfo("ENDVARIANCE", (oldEndVar + newEnd * newEnd) + "");
				
				// Update SVLEN
				long oldLength = consensus[groupNumber].getLength();
				long newLength = entry.getLength();
				consensus[groupNumber].setInfo("SVLEN", (oldLength + newLength) + "");
				
				String varId = entry.getId();
				varId = varId.substring(varId.indexOf('_') + 1);
				idLists[groupNumber].append("," + varId);
				
				if(Settings.INPUTS_MERGED)
				{
					// Update the extended support vector
					String suppVecExt = entry.getInfo("SUPP_VEC_EXT");
					if(suppVecExt.length() == 0)
					{
						suppVecExt = entry.getInfo("SUPP_VEC");
					}
					if(suppVecExt.length() == 0)
					{
						suppVecExt = "1";
					}
					
					// Add zeroes for any absent sample before this
					for(int i = sample-1; i>=0; i--)
					{
						if(supportVectors[groupNumber].charAt(i) == '1')
						{
							break;
						}
						int previouslyMergedCount = VariantInput.previouslyMergedSamples.get(i);
						for(int j = 0; j<previouslyMergedCount; j++)
						{
							suppVecExt = "0" + suppVecExt;
						}
					}
					suppVecExt = consensus[groupNumber].getInfo("SUPP_VEC_EXT") + suppVecExt;
					consensus[groupNumber].setInfo("SUPP_VEC_EXT", suppVecExt);
					
					// Update the extended ID List
					String idListExt = entry.getInfo("IDLIST_EXT");
					if(idListExt.length() == 0)
					{
						idListExt = entry.getInfo("IDLIST");
					}
					if(idListExt.length() == 0)
					{
						idListExt = entry.oldId;
					}
					
					String oldIdListExt = consensus[groupNumber].getInfo("IDLIST_EXT");
					
					consensus[groupNumber].setInfo("IDLIST_EXT", oldIdListExt + "," + idListExt);
					
					// Update specific marker
					if(consensus[groupNumber].hasInfoField("IS_SPECIFIC") && entry.hasInfoField("IS_SPECIFIC"))
					{
						if(consensus[groupNumber].getInfo("IS_SPECIFIC").equals("0") && entry.getInfo("IS_SPECIFIC").equals("1"))
						{
							consensus[groupNumber].setInfo("IS_SPECIFIC", "1");
						}
					}
					
					// Update extended support count
					int extendedSupport = 1;
					int oldExtendedSupport = Integer.parseInt(consensus[groupNumber].getInfo("SUPP_EXT"));
					if(entry.hasInfoField("SUPP_EXT"))
					{
						extendedSupport = Integer.parseInt(entry.getInfo("SUPP_EXT"));
					}
					else if(entry.hasInfoField("SUPP"))
					{
						extendedSupport = Integer.parseInt(entry.getInfo("SUPP"));
					}
					consensus[groupNumber].setInfo("SUPP_EXT", (extendedSupport + oldExtendedSupport) + "");
				}
			}
			
			
			used[groupNumber]++;
			
			// If this group is done, divide out any averages (e.g., position) as necessary
			if(used[groupNumber] == sizes[groupNumber])
			{
				if((!Settings.USE_STRAND))
				{
					consensus[groupNumber].setInfo("STRANDS", "??");
				}
				if(!Settings.USE_TYPE)
				{
					consensus[groupNumber].setInfo("SVTYPE", "???");
				}
				
				// Divide start and end by number of merged variants and round
				int groupSize = sizes[groupNumber];
				long totalStart = consensus[groupNumber].getPos();
				long totalEnd = consensus[groupNumber].getEnd();
				long totalSquaredStart = Long.parseLong(consensus[groupNumber].getInfo("STARTVARIANCE"));
				long totalSquaredEnd = Long.parseLong(consensus[groupNumber].getInfo("ENDVARIANCE"));
				
				// Compute start and end variances from the INFO fields
				double expectedSquaredStart = totalSquaredStart * 1.0 / groupSize;
				double expectedStartSquared = totalStart * totalStart * 1.0 / groupSize / groupSize;
				double varStart = expectedSquaredStart - expectedStartSquared;
				double expectedSquaredEnd = totalSquaredEnd * 1.0 / groupSize;
				double expectedEndSquared = totalEnd * totalEnd * 1.0 / groupSize / groupSize;
				double varEnd = expectedSquaredEnd - expectedEndSquared;
				String varStartString = String.format("%.6f", varStart);
				String varEndString = String.format("%.6f", varEnd);
				
				// Update the start and end mean/variance values
				consensus[groupNumber].setPos((long)((totalStart * 1.0 + .5) / groupSize)); 
				consensus[groupNumber].setInfo("END", (long)((totalEnd * 1.0 + .5) / groupSize) + "");
				consensus[groupNumber].setInfo("STARTVARIANCE", varStartString);
				consensus[groupNumber].setInfo("ENDVARIANCE", varEndString);
				
				// Update the average SV length
				long totalLength = consensus[groupNumber].getLength();
				consensus[groupNumber].setInfo("SVLEN", (long)((totalLength * 1.0 + .5) / groupSize) + "");
				
				// Fill the support-related fields
				consensus[groupNumber].setInfo("SUPP_VEC", supportVectors[groupNumber]);
				consensus[groupNumber].setInfo("SUPP", supportCounts[groupNumber]+"");
				consensus[groupNumber].setInfo("SVMETHOD", "JASMINE");
				consensus[groupNumber].setInfo("IDLIST", idLists[groupNumber].toString());
				
				if(Settings.KEEP_FIRST_POS)
				{
					consensus[groupNumber].setPos(firstStart.get(groupNumber));
					consensus[groupNumber].setInfo("END", firstEnd.get(groupNumber) + "");
					consensus[groupNumber].setInfo("SVLEN", firstLength.get(groupNumber) + "");
				}
				
				if(Settings.INPUTS_MERGED)
				{
					// Add zeroes to SUPP_VEC_EXT as needed
					for(int i = sample+1; i < supportVectors[groupNumber].length(); i++)
					{
						int previouslyMergedCount = VariantInput.previouslyMergedSamples.get(i);
						for(int j = 0; j<previouslyMergedCount; j++)
						{
							consensus[groupNumber].setInfo("SUPP_VEC_EXT", consensus[groupNumber].getInfo("SUPP_VEC_EXT") + "0");
						}
					}
				}
				// Remove the sample number from the variant ID (copied over from the first sample which is a part of this merged set)
				if(!Settings.CHANGE_VAR_IDS)
				{
					consensus[groupNumber].setId(consensus[groupNumber].oldId);
				}
				
				if(supportCounts[groupNumber] >= Settings.MIN_SUPPORT)
				{
					out.println(consensus[groupNumber]);
				}
				consensus[groupNumber] = null;
			}
		}
		
		
	}
}
