/*
 * Script for visualizing all variants in a merged VCF file
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;

public class IgvScreenshotMaker {
	
	static String vcfFn = "";
	static String bamFilelist = "";
	static String vcfFilelist = "";
	static String genomeFn = "";
	
	static String outPrefix = "";
	
	static boolean SQUISH = false;
	static boolean SVG = false;
	static boolean PRECISE = false;
	
	static HashMap<String, String> infoFilters;
	static HashSet<String> grepFilters;
		
	static void parseArgs(String[] args)
	{
		infoFilters = new HashMap<String, String>();
		grepFilters = new HashSet<String>();
		
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf('=');
			if(equalsIdx == -1)
			{
				if(arg.toLowerCase().endsWith("squish"))
				{
					SQUISH = true;
				}
				else if(arg.toLowerCase().endsWith("svg"))
				{
					SVG = true;
				}
				else if(arg.toLowerCase().endsWith("normalize_chr_names"))
				{
					Settings.DEFAULT_CHR_NORM = true;
				}
				else if(arg.toLowerCase().endsWith("precise"))
				{
					PRECISE = true;
				}
			}
			else
			{
				String key = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("vcf_file"))
				{
					vcfFn = val;
				}
				else if(key.equalsIgnoreCase("genome_file"))
				{
					genomeFn = val;
				}
				else if(key.equalsIgnoreCase("bam_filelist"))
				{
					bamFilelist = val;
				}
				else if(key.equalsIgnoreCase("vcf_filelist"))
				{
					vcfFilelist = val;
				}
				else if(key.equalsIgnoreCase("out_prefix"))
				{
					outPrefix = val;
				}
				else if(key.equalsIgnoreCase("info_filter"))
				{
					String[] tokens = val.split(",");
					infoFilters.put(tokens[0], tokens[1]);
				}
				else if(key.equalsIgnoreCase("grep_filter"))
				{
					grepFilters.add(val);
				}
			}
		}
		
		if(vcfFn.length() == 0 || genomeFn.length() == 0 || bamFilelist.length() == 0 || outPrefix.length() == 0)
		{
			usage();
			System.exit(0);
		}
	}
	
	/*
	 * Print the usage menu
	 */
	static void usage()
	{
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  vcf_file      (String) - the VCF file with merged SVs");
		System.out.println("  genome_file   (String) - the FASTA file with the reference genome");
		System.out.println("  bam_filelist  (String) - a comma-separated list of BAM files");
		System.out.println("  out_prefix    (String) - the prefix of the output directory and filenames");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  info_filter=KEY,VALUE  - filter by an INFO field value (multiple allowed) e.g., info_filter=SUPP_VEC,101");
		System.out.println("  grep_filter=QUERY      - filter to only lines containing a given QUERY");
		System.out.println("  vcf_filelist  (String) - the txt file with a list of merged VCFs");
		System.out.println("  --precise              - require variant to contain \"PRECISE\" as an INFO field");
		System.out.println("  --squish               - squishes tracks to fit more reads");
		System.out.println("  --svg                  - save as an SVG instead of a PNG");
		System.out.println("  --normalize_chr_names  - normalize the VCF chromosome name to strip \"chr\"");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		Settings.CHR_NAME_MAP = new ChrNameNormalization();
		
		parseArgs(args);
		
		Path currentRelativePath = Paths.get("");
		String outDir = currentRelativePath.toAbsolutePath().toString() + "/" + outPrefix;
		File outDirFile = new File(outDir);
		if(outDirFile.isDirectory())
		{
			final File[] files = outDirFile.listFiles();
			for (File f: files) f.delete();
			outDirFile.delete();
		}
		outDirFile.mkdir();
		String ofn = outDir + "/" + outPrefix + ".bat";
		
		PrintWriter out = new PrintWriter(new File(ofn));
		
		out.println("new");
		out.println("genome " + (genomeFn.startsWith("/") ? 
				genomeFn : (currentRelativePath.toAbsolutePath().toString() + "/" + genomeFn)));
		ArrayList<String> bamFiles = PipelineManager.getFilesFromList(bamFilelist);
		ArrayList<String> vcfFiles = new ArrayList<String>();
		if(vcfFilelist.length() > 0)
		{
			vcfFiles = PipelineManager.getFilesFromList(vcfFilelist);
		}
		for(int i = 0; i<bamFiles.size(); i++)
		{
			String bamFn = bamFiles.get(i);
			out.println("load " + (bamFn.startsWith("/") ? 
					bamFn : (currentRelativePath.toAbsolutePath().toString() + "/" + bamFn)));
			if(vcfFiles.size() > 0)
			{
				String fn = currentRelativePath.toAbsolutePath().toString() + "/" + StringUtils.fileBaseName(bamFn);
				fn = fn.substring(0, fn.length() - 4) + ".bed";
				out.println("load " + fn);
				
				PrintWriter curOut = new PrintWriter(new File(fn));
				Scanner curInput = new Scanner(new FileInputStream(new File(vcfFiles.get(i))));
				while(curInput.hasNext())
				{
					String line = curInput.nextLine();
					if(line.length() == 0 || line.startsWith("#"))
					{
						continue;
					}
					VcfEntry entry = VcfEntry.fromLine(line);
					String chr = entry.getChromosome();
					int start = (int)entry.getPos();
					int end = (int)entry.getEnd();
					String id = entry.getId();
					String type = entry.getNormalizedType();
					if(type.equalsIgnoreCase("TRA"))
					{
						String chr2 = entry.getChr2();
						curOut.printf("%s\t%d\t%d\t%s_%s\n", chr, start, start+1, id, type);
						curOut.printf("%s\t%d\t%d\t%s_%s\n", chr2, end, end+1, id, type);
					}
					else
					{
						if(end - start <= 100000)
						{
							curOut.printf("%s\t%d\t%d\t%s_%s\n", chr, start, end+1, id, type);
						}
					}
				}
				curInput.close();
				curOut.close();
			}
			
		}
		out.println("snapshotDirectory " + outDir);
		
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0 || line.startsWith("#"))
			{
				continue;
			}
			VcfEntry entry = new VcfEntry(line);
			
			// Check that the entry passes grep and INFO filters
			boolean passesFilters = true;
			
			for(String s : grepFilters)
			{
				if(!line.contains(s))
				{
					passesFilters = false;
				}
			}
			for(String s : infoFilters.keySet())
			{
				if(!entry.hasInfoField(s) || !entry.getInfo(s).equals(infoFilters.get(s)))
				{
					passesFilters = false;
				}
			}
			
			if(PRECISE && !entry.tabTokens[7].startsWith("PRECISE;") && !entry.tabTokens[7].contains(";PRECISE;"))
			{
				passesFilters = false;
			}
			
			if(!passesFilters)
			{
				continue;
			}
			
			long start = entry.getPos() - 100;
			long end = entry.getEnd() + 100;
			if(entry.getNormalizedType().equals("INS"))
			{
				end = start + entry.getLength() + 100;
			}
			
			if(end > start + 100000)
			{
				continue;
			}
			
			String chr = entry.getChromosome();
			
			out.println("goto " + chr + ":" + start + "-" + end);
			out.println("sort position");
			if(SQUISH)
			{
				for(String bamFile : bamFiles)
				{
					out.println("squish " + bamFile);
				}
			}
			else
			{
				for(String bamFile : bamFiles)
				{
					out.println("collapse " + bamFile);
				}
			}
			out.println("snapshot " + entry.getId() + ".png");
		}
		
		out.println("exit");
		
		input.close();
		out.close();
		
		out.close();
	}
}
