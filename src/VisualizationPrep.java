/*
 * Given the results of both THRIVER and SURVIVOR, extract out shared and different merges,
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
	@SuppressWarnings("unchecked")
	public static void main(String[] args) throws Exception
	{
		//String fileList = "/home/mkirsche/eichler/filelist.txt";
		String fileList = "filelist.txt";
		//String survivorOutput = "/home/mkirsche/eichler/survmerged.vcf";
		//String thriverOutput = "/home/mkirsche/eichler/merged.vcf";
		String thriverOutput = "out.vcf";
		String survivorOutput = "outsurv.vcf";
		Scanner fileNameReader = new Scanner(new FileInputStream(new File(fileList)));
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
		
		String outFile = thriverOutput + ".graph";
		
		int[] ys = new int[vcfs.length];
		for(int i = 0; i<ys.length; i++)
		{
			ys[i] = i;
		}
		
		PrintWriter out = new PrintWriter(new File(outFile));
		
		// For each variant, keep track of its position and color (type)
		ArrayList<Integer>[] positions = new ArrayList[vcfs.length];
		ArrayList<Integer>[] colors = new ArrayList[vcfs.length];
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
				if(!entry.getChromosome().equals("1")) continue;
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
		
		// Now we have to get line segments, so get merged sets from both SURVIVOR and THRIVER.
		TreeSet<Merge> myEdges = getJoinedPairs(thriverOutput, false);
		TreeSet<Merge> survEdges = getJoinedPairs(survivorOutput, true);
		
		// Store the union of the merge-sets so we get every line segment
		TreeSet<Merge> union = new TreeSet<Merge>();
		for(Merge s : myEdges) union.add(s);
		for(Merge s : survEdges) union.add(s);
		
		// For each line segment, color it based on which software it came from (possibly both)
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
			if(myEdges.contains(edge)) color |= 1;
			if(survEdges.contains(edge)) color |= 2;
			
			// If the pair was only merged by one software, print out information about it
			if(color == 1 || color == 2)
			{
				VcfEntry first = idToEntry[samples[0]].get(ids[0]);
				VcfEntry second = idToEntry[samples[1]].get(ids[1]);
				System.out.println("Merge unique to " + (color == 1 ? "thriver" : "survivor"));
				System.out.println("  " + ids[0] + " " + first.getType() + " " + first.getStrand() + " at " + first.getPos() + " (length " + first.getLength() + ")");
				System.out.println("  " + ids[1] + " " + second.getType() + " " + second.getStrand() + " at " + second.getPos() + " (length " + second.getLength() + ")");
				System.out.println("  " + edge.line);
				System.out.println("  Samples: " + edge.sample1 + " " + edge.sample2);
				Variant a = VariantInput.fromVcfEntry(first, 0), b = VariantInput.fromVcfEntry(second, 0);
				System.out.println("  Distance according to THRIVER: " + a.distance(b));
			}
			
			// Print the line segment
			out.println(curPositions[0]+" "+ys[edge.sample1]+" "+curPositions[1]+" "+ys[edge.sample2]+" "+color);
		}
		out.close();
		
		
	}
	
	/*
	 * For a given merged VCF file, get the list of all pairs of variants which were joined
	 * For now, assumes only 2 samples, and the survivor flag is true if SURVIVOR was used
	 * and false if THRIVER was used instead.
	 */
	static TreeSet<Merge> getJoinedPairs(String fn, boolean survivor) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		TreeSet<Merge> res = new TreeSet<Merge>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#")) continue;
			VcfEntry entry = new VcfEntry(line);
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
				//String[] ids = new String[entry.tabTokens.length - 9];
				ArrayList<String> ids = new ArrayList<String>();
				for(int i = 9; i<entry.tabTokens.length; i++)
				{
					String val = entry.tabTokens[i].split(":")[7];
					if(!val.equalsIgnoreCase("nan")) ids.add(val);
					//ids.add();
					//ids[i] = entry.tabTokens[i+9].split(":")[7];
				}
				for(int i = 0; i<ids.size()-1 && i < samples.size()-1; i++)
				{
					res.add(new Merge(ids.get(i), ids.get(i+1), samples.get(i), samples.get(i+1), line));
				}
				
			}
			else
			{
				String[] ids = entry.getInfo("IDLIST").split(",");
				for(int i = 0; i<ids.length-1 && i < samples.size()-1; i++)
				{
					res.add(new Merge(ids[i], ids[i+1], samples.get(i), samples.get(i+1), line));
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
