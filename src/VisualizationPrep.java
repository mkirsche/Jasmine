/*
 * Given the results of both Jasmine and SURVIVOR, extract out shared and different merges,
 * producing a list of points and line segments which can be plotted to visualize the results.
 * 
 * For now, this only works on datasets with two samples (i.e., two VCFs input to the merging software).
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeSet;
public class VisualizationPrep {
	
	// Empty string if whole genome or chromosome name for plotting that chromosome
	static String chrToPlot = "1";
	
	 // True iff we did one of the two with the samples in reverse order
	static boolean secondRev = true;
	
	// Whether or not each file was produced by SURVIVOR and needs to be parsed differently
	static boolean firstSurvivor = false;
	static boolean secondSurvivor = false;
	
	// Whether or not to print merges unique to one output file
	static boolean printUnique = false;
	
	static int sampleCount = 0;
	
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception
	{
		// File containing a list of VCF files
		String fileList = "/home/mkirsche/eichler/filelist.txt";
		
		// The resulting merged VCF files from both Jasmine and SURVIVOR
		String firstOutput = "/home/mkirsche/eichler/merged.vcf";
		String secondOutput = "/home/mkirsche/eichler/revmerged.vcf";
		Scanner fileNameReader = new Scanner(new FileInputStream(new File(fileList)));
		
		// Get the list of VCF files
		ArrayList<String> vcfsAsList = new ArrayList<String>();
		while(fileNameReader.hasNext())
		{
			String line = fileNameReader.nextLine();
			if(line.length() == 0) continue;
			vcfsAsList.add(line);
		}
		fileNameReader.close();
		String[] vcfs = new String[vcfsAsList.size()];
		for(int i = 0; i<vcfs.length; i++)
		{
			vcfs[i] = vcfsAsList.get(i);
		}
		
		sampleCount = vcfs.length;
		
		String outFile = secondOutput + ".graph";
		
		// The y-coordinate of each point (variant) will be the sample it came from
		int[] ys = new int[vcfs.length];
		for(int i = 0; i<ys.length; i++)
		{
			ys[i] = i;
		}
		
		PrintWriter out = new PrintWriter(new File(outFile));
		
		// For each variant, keep track of its position and color (type)
		ArrayList<Integer>[] positions = new ArrayList[vcfs.length];
		ArrayList<Integer>[] colors = new ArrayList[vcfs.length];
		int[] colorCounts = new int[4];
		HashMap<String, VcfEntry>[] idToEntry = new HashMap[vcfs.length];
		
		// Hard-code the colors of the common variant types, so we know which color is which in downstream plotting
		// There may be other colors in the case of other types.
		HashMap<String, Integer> typeToInt = new HashMap<String, Integer>();
		typeToInt.put("INS", 0);
		typeToInt.put("DEL", 1);
		typeToInt.put("DUP", 2);
		typeToInt.put("INV", 3);
		
		// Iterate over input VCF files to record entries
		for(int i = 0; i<vcfs.length; i++)
		{
			positions[i] = new ArrayList<Integer>();
			colors[i] = new ArrayList<Integer>();
			idToEntry[i] = new HashMap<String, VcfEntry>();
			Scanner input = new Scanner(new FileInputStream(new File(vcfs[i])));
			
			// Read entries one at a time and store information about them
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.length() == 0 || line.startsWith("#"))
				{
					continue;
				}
				VcfEntry entry = new VcfEntry(line);
				if(chrToPlot.length() > 0 && !entry.getChromosome().equals(chrToPlot)) continue;
				
				// Below is an example of how to restrict the plot to certain positions
				//if(entry.getPos() > 10000000) continue;
				
				int pos = (int)entry.getPos();
				positions[i].add(pos);
				idToEntry[i].put(entry.getId(), entry);
				if(!typeToInt.containsKey(entry.getType()))
				{
					typeToInt.put(entry.getType(), typeToInt.size());
				}
				colors[i].add(typeToInt.get(entry.getType()));
			}
			input.close();
		}
		
		// We now have the variant points, so output those for plotting
		for(int i = 0; i<vcfs.length; i++)
		{
			for(int j = 0; j<positions[i].size(); j++)
			{
				int x = positions[i].get(j);
				out.println(x+" "+ys[i]+" "+colors[i].get(j));
			}
		}
		
		// Now we have to get line segments, so get merged sets from both SURVIVOR and Jasmine.
		TreeSet<Merge> firstEdges = getJoinedPairs(firstOutput, firstSurvivor, false);
		System.out.println("Merges in first output: " + firstEdges.size());
		TreeSet<Merge> secondEdges = getJoinedPairs(secondOutput, secondSurvivor, secondRev);
		System.out.println("Merges in second output: " + secondEdges.size());
		
		// Store the union of the merge-sets so we get every line segment
		TreeSet<Merge> union = new TreeSet<Merge>();
		for(Merge s : secondEdges) union.add(s);
		for(Merge s : firstEdges) union.add(s);
		
		// For each line segment, color it based on which output it came from (possibly both)
		for(Merge edge : union)
		{
			String[] ids = new String[] {edge.id1, edge.id2};
			int[] samples = new int[] {edge.sample1, edge.sample2};
			boolean okay = true;
			int[] curPositions = new int[2];
			for(int i = 0; i<2; i++)
			{
				if(idToEntry[samples[i]].containsKey(ids[i]))
				{
					curPositions[i] = (int)idToEntry[samples[i]].get(ids[i]).getPos();
				}
				else
				{
					okay = false;
					break;
				}
			}
			if(!okay)
			{
				continue;
			}
			int color = 0;
			if(secondEdges.contains(edge)) color |= 2;
			if(firstEdges.contains(edge)) color |= 1;
			
			String firstSoftware = firstSurvivor ? "survivor" : "Jasmine";
			String secondSoftware = secondSurvivor ? "survivor" : "Jasmine";
			if(secondRev) secondSoftware += "rev";
			
			colorCounts[color]++;
			
			// If the pair was only merged by one software, print out information about it
			if(color == 1 || color == 2)
			{
				VcfEntry first = idToEntry[samples[0]].get(ids[0]);
				VcfEntry second = idToEntry[samples[1]].get(ids[1]);
				
				if(printUnique)
				{
					System.out.println("Merge unique to " + (color == 1 ? firstSoftware : secondSoftware));
					System.out.println("  " + ids[0] + " " + first.getType() + " " + first.getStrand() + " at " + first.getPos() + " (length " + first.getLength() + ")");
					System.out.println("  " + ids[1] + " " + second.getType() + " " + second.getStrand() + " at " + second.getPos() + " (length " + second.getLength() + ")");
					System.out.println("  " + edge.line);
					System.out.println("  Samples: " + edge.sample1 + " " + edge.sample2);
					Variant a = VariantInput.fromVcfEntry(first, 0), b = VariantInput.fromVcfEntry(second, 0);
					System.out.println("  Distance according to Jasmine: " + a.distance(b));
				}
			}
			
			// Print the line segment
			out.println(curPositions[0]+" "+ys[edge.sample1]+" "+curPositions[1]+" "+ys[edge.sample2]+" "+color);
		}
		System.out.println("First output unique merges: " + colorCounts[1]);
		System.out.println("Second output unique merges: " + colorCounts[2]);
		System.out.println("Shared merges: " + colorCounts[3]);
		out.close();
		
		
	}
	
	/*
	 * For a given merged VCF file, get the list of all pairs of variants which were joined
	 * For now, assumes only 2 samples, and the survivor flag is true if SURVIVOR was used
	 * and false if Jasmine was used instead.
	 */
	static TreeSet<Merge> getJoinedPairs(String fn, boolean survivor, boolean rev) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		TreeSet<Merge> res = new TreeSet<Merge>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#")) continue;
			VcfEntry entry = new VcfEntry(line);
			if(chrToPlot.length() > 0 && !entry.getChromosome().equals(chrToPlot)) continue;
			String supportVector = entry.getInfo("SUPP_VEC");
			ArrayList<Integer> samples = new ArrayList<Integer>();
			for(int i = 0; i<supportVector.length(); i++)
			{
				if(supportVector.charAt(i) == '1')
				{
					samples.add(i);
				}
			}
			if(survivor)
			{
				if(entry.tabTokens.length < 11) continue;
				ArrayList<String> ids = new ArrayList<String>();
				for(int i = 9; i<entry.tabTokens.length; i++)
				{
					String val = entry.tabTokens[i].split(":")[7];
					if(!val.equalsIgnoreCase("nan")) ids.add(val);
				}
				for(int i = 0; i<ids.size()-1 && i < samples.size()-1; i++)
				{
					if(rev) res.add(new Merge(ids.get(i+1), ids.get(i), sampleCount - 1 - samples.get(i+1), sampleCount - 1 - samples.get(i), line));
					//if(rev) res.add(new Merge(ids.get(i), ids.get(i+1), samples.get(samples.size()-1-i), samples.get(samples.size()-1-(i+1)), line));
					else res.add(new Merge(ids.get(i), ids.get(i+1), samples.get(i), samples.get(i+1), line));
				}
				
			}
			else
			{
				String[] ids = entry.getInfo("IDLIST").split(",");
				for(int i = 0; i<ids.length-1 && i < samples.size()-1; i++)
				{
					if(rev) res.add(new Merge(ids[i+1], ids[i], sampleCount - 1 - samples.get(i+1), sampleCount - 1 - samples.get(i), line));
					else res.add(new Merge(ids[i], ids[i+1], samples.get(i), samples.get(i+1), line));
				}
			}
			
		}
		input.close();
		return res;
	}
	
	/*
	 * A merge between two variants in different samples
	 */
	static class Merge implements Comparable<Merge>
	{
		String id1, id2;
		int sample1, sample2;
		String line;
		Merge(String ii1, String ii2, int ss1, int ss2, String ll)
		{
			line = ll;
			id1 = ii1;
			id2 = ii2;
			sample1 = ss1;
			sample2 = ss2;
		}
		@Override
		public int compareTo(Merge o) {
			if(sample1 != o.sample1)
			{
				return sample1 - o.sample1;
			}
			if(sample2 != o.sample2)
			{
				return sample2 - o.sample2;
			}
			if(!id1.equals(o.id1))
			{
				return id1.compareTo(o.id1);
			}
			return id2.compareTo(o.id2);
		}
	}
}
