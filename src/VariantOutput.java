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
	
	void addGraph(String graphID, ArrayList<Variant>[] graph, int sampleCount)
	{
		groups.put(graphID,  new VariantGraph(graph, sampleCount));
	}
	
	// Given a list of VCF files and merging results, output an updated VCF file
	// TODO use minSupport
	public void writeMergedVariants(String fileList, String outFile, int minSupport) throws Exception
	{
		Scanner listInput = new Scanner(new FileInputStream(new File(fileList)));
		PrintWriter out = new PrintWriter(new File(outFile));
		int sample = 0;
		
		ArrayList<String> headerLinesToPrint = new ArrayList<String>();
		
		// Go through one VCF file at a time
		while(listInput.hasNext())
		{
			String filename = listInput.nextLine();
			if(filename.length() == 0)
			{
				continue;
			}
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
						headerLinesToPrint.add(line);
					}
					else
					{
						continue;
					}
				}
				else
				{
					// Update the consensus variant in the appropriate graph
					VcfEntry entry = new VcfEntry(line);
					String graphID = entry.getGraphID();
					groups.get(graphID).processVariant(entry, sample);
				}
			}
			input.close();
			sample++;
		}
		
		// Go through the graphs and add all the consensus variants to a big list
		ArrayList<VcfEntry> allEntries = new ArrayList<VcfEntry>();
		for(String graphID : groups.keySet())
		{
			VariantGraph graph = groups.get(graphID);
			for(int i = 0; i<graph.consensus.length; i++)
			{
				if(graph.consensus[i] == null)
				{
					continue;
				}
				if(graph.supportCounts[i] >= minSupport)
				{
					allEntries.add(graph.consensus[i]);
				}
			}
		}
		
		// Actually print the header and variants
		for(String headerLine : headerLinesToPrint)
		{
			out.println(headerLine);
		}
		for(VcfEntry entry : allEntries)
		{
			String oldId = entry.getId();
			entry.setId(oldId.substring(1 + oldId.indexOf('_')));
			out.println(entry);
		}
		
		out.close();
		listInput.close();
	}
	
	static class VariantGraph
	{
		int[] sizes; // Number of variants in each group
		int[] used; // How many variants in each group have been seen so far
		VcfEntry[] consensus; // The current consensus variant for each group
		HashMap<String, Integer> varToGroup; // For each variant ID, the group number it is in
		String[] supportVectors;
		int[] supportCounts;
		StringBuilder[] idLists; // The list of variant IDs in each merged variant
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
		
		// From a VCF line, update the appropriate consensus entry
		void processVariant(VcfEntry entry, int sample) throws Exception
		{
			// This should never happen, but if the variant ID is not in the graph ignore it
			String fullId = VariantInput.fromVcfEntry(entry, sample).id;
			if(!varToGroup.containsKey(fullId))
			{
				return;
			}
			
			int groupNumber = varToGroup.get(fullId);
			
			// If this is the first variant in the group, initialize the consensus entry
			if(used[groupNumber] == 0)
			{
				consensus[groupNumber] = entry;
				String varId = entry.getId();
				varId = varId.substring(varId.indexOf('_') + 1);
				idLists[groupNumber].append(varId);
			}
			
			// Otherwise, update the consensus to include info from this variant
			// For average, don't divide yet to avoid loss of precision
			else
			{
				consensus[groupNumber].setPos(consensus[groupNumber].getPos() + entry.getPos());
				String varId = entry.getId();
				varId = varId.substring(varId.indexOf('_') + 1);
				idLists[groupNumber].append("," + varId);
			}
			
			used[groupNumber]++;
			
			// If this group is done, divide out any averages (e.g., position) as necessary
			if(used[groupNumber] == sizes[groupNumber])
			{
				consensus[groupNumber].setPos((long)(.5 + 1.0 * consensus[groupNumber].getPos() / sizes[groupNumber]));
				consensus[groupNumber].setInfo("SUPP_VEC", supportVectors[groupNumber]);
				consensus[groupNumber].setInfo("SUPP", supportCounts[groupNumber]+"");
				consensus[groupNumber].setInfo("SVMETHOD", "THRIVER");
				consensus[groupNumber].setInfo("IDLIST", idLists[groupNumber].toString());
				String varId = entry.getId();
				varId = varId.substring(varId.indexOf('_') + 1);
				consensus[groupNumber].setId(varId);
			}
		}
		
		
	}
}
