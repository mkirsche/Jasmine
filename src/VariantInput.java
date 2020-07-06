/*
 * Methods for reading VCF entries from VCF files and dividing
 * the entries into separate groups by graph ID
 */

import java.io.File;
import java.io.FileInputStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeMap;

public class VariantInput {
	
	// How many samples were merged to produce each input file
	static HashMap<Integer, Integer> previouslyMergedSamples = new HashMap<Integer, Integer>();
	
	/*
	 * Count the number of VCF files in a list
	 */
	public static int countFiles(String fileList) throws Exception
	{
		return PipelineManager.getFilesFromList(fileList).size();
	}
	
	/*
	 * Get a list of all variants from a group of files, binning them by graphID
	 */
	@SuppressWarnings("unchecked")
	public static TreeMap<String, ArrayList<Variant>> readAllFiles(String fileList) throws Exception
	{
		ArrayList<String> fileNames = PipelineManager.getFilesFromList(fileList);
		
		TreeMap<String, ArrayList<Variant>>[] variantsPerFile = new TreeMap[fileNames.size()];
		for(int i = 0; i<fileNames.size(); i++)
		{
			variantsPerFile[i] = getSingleList(fileNames.get(i), i);
		}
		TreeMap<String, ArrayList<Variant>> res = new TreeMap<String, ArrayList<Variant>>();
		for(int i = 0; i<fileNames.size(); i++)
		{
			for(String s : variantsPerFile[i].keySet())
			{
				if(!res.containsKey(s))
				{
					res.put(s, new ArrayList<Variant>());
				}
				for(Variant v : variantsPerFile[i].get(s))
				{
					res.get(s).add(v);
				}
			}
		}
		return res;
	}
	
	/*
	 * Get a list of variants binned by graphID for a single VCF file
	 */
	private static TreeMap<String, ArrayList<Variant>> getSingleList(String filename, int sample) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(filename)));
		ArrayList<Variant> allVariants = new ArrayList<Variant>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			VcfEntry entry = VcfEntry.fromLine(line);
			if(!previouslyMergedSamples.containsKey(sample))
			{
				if(entry.getInfo("SUPP_VEC_EXT").length() > 0)
				{
					previouslyMergedSamples.put(sample, entry.getInfo("SUPP_VEC_EXT").length());
				}
				else if(entry.getInfo("SUPP_VEC").length() > 0)
				{
					previouslyMergedSamples.put(sample, entry.getInfo("SUPP_VEC").length());
				}
				else
				{
					previouslyMergedSamples.put(sample, 1);
				}
			}
			allVariants.add(fromVcfEntry(entry, sample));
			
		}
		
		System.out.println(filename + " has " + allVariants.size() + " variants");
		input.close();
		
		return divideIntoGraphs(allVariants);
	}
	
	/*
	 * Take a list of variants and bin them by graphID
	 */
	private static TreeMap<String, ArrayList<Variant>> divideIntoGraphs(ArrayList<Variant> data)
	{
		TreeMap<String, ArrayList<Variant>> groups = new TreeMap<String, ArrayList<Variant>>();
		for(Variant v : data)
		{
			String graphID = v.graphID;
			if(!groups.containsKey(graphID))
			{
				groups.put(graphID, new ArrayList<Variant>());
			}
			groups.get(graphID).add(v);
		}
		return groups;
	}
	
	/*
	 * From a line of a VCF file, extract the information needed for merging
	 * and return it as a Variant object
	 */
	public static Variant fromVcfEntry(VcfEntry entry, int sample) throws Exception
	{
		double start = entry.getFirstCoord();
		double end = entry.getSecondCoord();
		
		entry.setId(sample + "_" + entry.getId());
		
		String id = entry.getGraphID();
		
		String seq = null;
		if(entry.getType().equals("INS"))
		{
			String entrySeq = entry.getSeq();
			if(entrySeq.length() > 0)
			{
				seq = entrySeq;
			}
		}
		
		// Default distance threshold model is constant, so set to that first
		int maxDist = Settings.MAX_DIST;
		double minSeqId = Settings.MIN_SEQUENCE_SIMILARITY;
		
		// Then, check if there is a per-variant distance threshold
		String maxDistInfo = entry.getInfo("JASMINE_DIST");
		if(maxDistInfo.length() > 0)
		{
			maxDist = Integer.parseInt(maxDistInfo);
		}
		
		// Next, check if there is a length-based threshold
		else if(Settings.MAX_DIST_LINEAR > 0)
		{
			maxDist = (int)(Settings.MAX_DIST_LINEAR * Math.abs(entry.getLength()) + 0.5);
			if(Settings.MAX_DIST_SET)
			{
				maxDist = Math.min(maxDist, Settings.MAX_DIST);
			}
			if(Settings.MIN_DIST != -1)
			{
				maxDist = Math.max(maxDist, Settings.MIN_DIST);
			}
		}
		
		// Check for per-variant sequence ID thresholds
		String minIdInfo = entry.getInfo("JASMINE_ID");
		if(minIdInfo.length() > 0)
		{
			minSeqId = Double.parseDouble(minIdInfo);
		}
		
		Variant res = new Variant(sample, entry.getId(), start, end, id, seq, maxDist, minSeqId);
		res.hash = Variant.hash(entry.tabTokens[7]);
		return res;
	}
}
