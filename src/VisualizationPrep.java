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
		String vcf1 = "/home/mkirsche/crossstitch/ENC0002.spes.rck.vcf";
		String vcf2 = "/home/mkirsche/crossstitch/ENC0003.spes.rck.vcf";
		String thriverOutput = "out.vcf";
		String survivorOutput = "outsurv.vcf";
		
		String outFile = thriverOutput + ".graph";
		
		String[] vcfs = new String[] {vcf1, vcf2};
		int[] ys = new int[] {0, 1}; // The y-coordinate where we want to plot each sample's variants
		
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
		TreeSet<String> myEdges = getJoinedPairs(thriverOutput, false);
		TreeSet<String> survEdges = getJoinedPairs(survivorOutput, true);
		
		// Store the union of the merge-sets so we get every line segment
		TreeSet<String> union = new TreeSet<String>();
		for(String s : myEdges) union.add(s);
		for(String s : survEdges) union.add(s);
		
		// For each line segment, color it based on which software it came from (possibly both)
		for(String edge : union)
		{
			String[] ids = edge.split(",");
			if(ids.length != 2)
			{
				continue;
			}
			boolean okay = true;
			int[] curPositions = new int[2];
			for(int i = 0; i<2; i++)
			{
				if(idToEntry[i].containsKey(ids[i]))
				{
					curPositions[i] = (int)idToEntry[i].get(ids[i]).getPos();
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
				VcfEntry first = idToEntry[0].get(ids[0]);
				VcfEntry second = idToEntry[1].get(ids[1]);
				System.out.println("Merge unique to " + (color == 1 ? "thriver" : "survivor"));
				System.out.println("  " + ids[0] + " " + first.getType() + " " + first.getStrand() + " at " + first.getPos() + " (length " + first.getLength() + ")");
				System.out.println("  " + ids[1] + " " + second.getType() + " " + second.getStrand() + " at " + second.getPos() + " (length " + second.getLength() + ")");
				Variant a = VariantInput.fromVcfEntry(first, 0), b = VariantInput.fromVcfEntry(second, 0);
				System.out.println("  Distance according to THRIVER: " + a.distance(b));
			}
			
			// Print the line segment
			out.println(curPositions[0]+" "+ys[0]+" "+curPositions[1]+" "+ys[1]+" "+color);
		}
		out.close();
		
		
	}
	
	/*
	 * For a given merged VCF file, get the list of all pairs of variants which were joined
	 * For now, assumes only 2 samples, and the survivor flag is true if SURVIVOR was used
	 * and false if THRIVER was used instead.
	 */
	static TreeSet<String> getJoinedPairs(String fn, boolean survivor) throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(fn)));
		TreeSet<String> res = new TreeSet<String>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#")) continue;
			VcfEntry entry = new VcfEntry(line);
			if(survivor)
			{
				if(entry.tabTokens.length < 11) continue;
				String[] ids = new String[] {
						entry.tabTokens[9].split(":")[7], 
						entry.tabTokens[10].split(":")[7]
				};
				res.add(ids[0] + "," + ids[1]);
			}
			else
			{
				String[] ids = entry.getInfo("IDLIST").split(",");
				res.add(ids[0]+ "," + ids[1]);
			}
			
		}
		input.close();
		return res;
	}
}
