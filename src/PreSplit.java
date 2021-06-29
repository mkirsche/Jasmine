import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class PreSplit
{
	static String fileList = "";
	static String outputDir = "";
	static int segmentLength = -1;
	static boolean transTogether = false;
	
	static void usage()
	{
		System.out.println();
		System.out.println("Usage: split_jasmine file_list output_dir [segment_length]");
		System.out.println("  Example: split_jasmine file_list=filelist.txt output_dir=/path/to/split_dir segment_length=10m");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  file_list       (String) - A txt file with a line-separated list of VCFs to be split");
		System.out.println("  output_dir      (String) - The directory to write the split files to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  segment_length  (int)    - length of segments to split chromosomes into (default whole-chromosome)");
		System.out.println("  --ignore_strand          - allow variants with different strands to be merged");
		System.out.println("  --ignore_type            - allow variants with different types to be merged");
		System.out.println("  --combine_translocations - keep all translocations together to reduce number of groups");
		System.out.println();
	}
	
	static void parseArgs(String[] args) throws Exception
	{
		for(int i = 0; i<args.length; i++)
		{
			if(args[i].indexOf('=') == -1)
			{
				if(args[i].endsWith("ignore_strand"))
				{
					Settings.USE_STRAND = false;
				}
				else if(args[i].endsWith("ignore_type"))
				{
					Settings.USE_TYPE = false;
				}
				else if(args[i].endsWith("combine_translocations"))
				{
					transTogether = true;
				}
			}
			else
			{
				int equalIdx = args[i].indexOf('=');
				String key = args[i].substring(0, equalIdx);
				while(key.length() > 0 && key.charAt(0) == '-')
				{
					key = key.substring(1);
				}
				String val = args[i].substring(1 + equalIdx);
				
				switch(key) 
				{
					case "segment_length":
						segmentLength = Settings.parseInt(val);
						break;
					case "file_list":
						fileList = val;
						break;
					case "output_dir":
						outputDir = val;
						break;
					default:
						break;
				}
			}
		}
		if(fileList.length() == 0 || outputDir.length() == 0)
		{
			usage();
			System.exit(0);
		}
		
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		String[] filelists = convertAll(fileList, outputDir, segmentLength);
		for(String s : filelists)
		{
			System.out.println(s);
		}
	}
	
	@SuppressWarnings("unchecked")
	static String[] convertAll(String fileList, String outDir, int segmentLength) throws Exception
	{
		if(!outDir.startsWith("/"))
		{
			Path currentRelativePath = Paths.get("");
			outDir = currentRelativePath.toAbsolutePath().toString() + "/" + outDir;
		}
		if(!new File(outDir).isDirectory())
		{
			new File(outDir).mkdir();
		}
		ArrayList<String> vcfFiles = PipelineManager.getFilesFromList(fileList);
		int n = vcfFiles.size();
		HashMap<String, String>[] splitMaps = new HashMap[n];
		HashSet<String> allKeys = new HashSet<String>();
		for(int i = 0; i<n; i++)
		{
			splitMaps[i] = convertFile(vcfFiles.get(i), outDir + "/splitSample" + i, segmentLength);
			allKeys.addAll(splitMaps[i].keySet());
		}
		String[][] splitList = new String[n][allKeys.size()];
		for(int i = 0; i<n; i++)
		{
			VcfHeader header = new VcfHeader();
			Scanner input = new Scanner(new FileInputStream(new File(vcfFiles.get(i))));
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
					continue;
				}
				break;
			}
			int idx = 0;
			for(String s : allKeys)
			{
				if(splitMaps[i].containsKey(s))
				{
					splitList[i][idx] = splitMaps[i].get(s);
				}
				else
				{
					String ofn = outDir + "/splitSample" + i + "_" + s + ".vcf";
					PrintWriter out = new PrintWriter(new File(ofn));
					header.print(out);
					out.close();
					splitList[i][idx] = ofn;
				}
				idx++;
			}
		}
		
		int idx = 0;
		String[] res = new String[allKeys.size()];
		for(String s : allKeys)
		{
			String ofn = outDir + "/" + "split_" + s + ".filelist.txt";
			PrintWriter out = new PrintWriter(new File(ofn));
			for(int i = 0; i<n; i++)
			{
				out.println(splitList[i][idx]);
			}
			out.close();
			res[idx] = ofn;
			idx++;
		}
		
		return res;
	}
	
	static HashMap<String, String> convertFile(String inputFile, String outputPrefix, int segmentLength) throws Exception
	{
		VcfHeader header = new VcfHeader();
		Scanner input = new Scanner(new FileInputStream(new File(inputFile)));
		HashMap<String, String> res = new HashMap<String, String>();
		HashMap<String, PrintWriter> writerMap = new HashMap<String, PrintWriter>();
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
				continue;
			}
			VcfEntry entry = VcfEntry.fromLine(line);
			String graphId = VariantInput.fromVcfEntry(entry, 0).graphID;
			if(segmentLength != -1 && !entry.getNormalizedType().equals("TRA"))
			{
				graphId = graphId + "_" + ((entry.getPos() / segmentLength) * segmentLength);
			}
			if(entry.getNormalizedType().equals("TRA") && transTogether)
			{
				graphId = "TRA";
			}
			if(!res.containsKey(graphId))
			{
				String ofn = outputPrefix + "_" + graphId + ".vcf";
				PrintWriter out = new PrintWriter(new File(ofn));
				header.print(out);
				res.put(graphId, ofn);
				writerMap.put(graphId, out);
			}
			PrintWriter out = writerMap.get(graphId);
			out.println(line);
		}
		for(String key : writerMap.keySet())
		{
			writerMap.get(key).close();
		}
		return res;
	}
}
