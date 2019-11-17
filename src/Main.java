/*
 * Main interface for Thriver.
 */

import java.util.ArrayList;
import java.util.TreeMap;

public class Main {
public static void main(String[] args) throws Exception
{
	Settings.parseArgs(args);
	
	// Get the variants and bin them into individual graphs
	TreeMap<String, ArrayList<Variant>> allVariants = VariantInput.readAllFiles(Settings.FILE_LIST);
	
	VariantOutput output = new VariantOutput();
	int totalMerged = 0;
	
	// Get the number of samples to know the length of the SUPP_VEC field
	int sampleCount = VariantInput.countFiles(Settings.FILE_LIST);
	
	// Merge one graph at a time
	for(String graphID : allVariants.keySet())
	{
		ArrayList<Variant> variantList = allVariants.get(graphID);
		VariantMerger vm = new VariantMerger(variantList);
		vm.runMerging();
		ArrayList<Variant>[] res = vm.getGroups();
		output.addGraph(graphID, res, sampleCount);
		int merges = 0;
		for(ArrayList<Variant> list : res)
		{
			if(list.size() > 1)
			{
				merges++;
			}
		}
		totalMerged += merges;
	}
	
	// Print the merged variants to a file if they have enough support
	output.writeMergedVariants(Settings.FILE_LIST, Settings.OUT_FILE, Settings.MIN_SUPPORT);
	System.out.println("Number of sets with multiple variants: " + totalMerged); 
	
}
}
