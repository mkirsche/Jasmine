import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;

/*
 * Utility for handling command line parameters and the usage message 
 */

public class Settings {
	
	static boolean USE_STRAND = true;
	static boolean USE_TYPE = true;
	static String FILE_LIST = "";
	static String OUT_FILE = "";
	static int MAX_DIST = 1000;
	static double MAX_DIST_LINEAR = 0.0;
	static int MIN_SUPPORT = 2;
	static double MIN_SEQUENCE_SIMILARITY = 0;
	
	static boolean CONVERT_DUPLICATIONS = false;
	static boolean MARK_SPECIFIC = false;
	static boolean RUN_IRIS = false;
	static String GENOME_FILE = "";
	static String BAM_FILE_LIST = "";
	static String IRIS_ARGS = "";
	
	static String OUT_DIR = "output";
	static int THREADS = 2;
	
	static int SPECIFIC_MIN_RCOUNT = 10;
	static int SPECIFIC_MIN_LENGTH = 30;
	
	static void usage()
	{
		System.out.println();
		System.out.println("Usage: java -cp src Main [args]");
		System.out.println("  Example: java -cp src Main file_list=filelist.txt out_file=out.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  file_list (String) - a file listing paths to all variant files to use (on separate lines)");
		System.out.println("  out_file  (String) - the name of the file to output the merged variants to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  max_dist        (int)    [1000] - the maximum distance variants can be apart when being merged");
		System.out.println("  max_dist_linear (float)  [0]    - make max_dist this proportion of the length of each variant (overrides max_dost)");
		System.out.println("  min_seq_id      (float)  [0]    - the minimum sequence identity for two insertions to be merged");
		System.out.println("  min_support     (int)    [2]    - the minimum number of callsets a variant must be in to be output");
		System.out.println("  threads         (int)    [2]    - the number of threads to use for merging the variants");
		System.out.println("  spec_reads      (int)    [10]   - the minimum number of reads a variant needs to be in the specific callset");
		System.out.println("  spec_len        (int)    [30]   - the minimum length a variant needs to be in the specific callset");
		System.out.println("  genome_file     (String) []     - the reference genome being used");
		System.out.println("  bam_list        (String) []     - a file listing paths to BAMs in the same order as the VCFs");
		System.out.println("  iris_args       (String) []     - a comma-separated list of optional arguments to pass to Iris");
		System.out.println("  --ignore_strand             - allow variants with different strands to be merged");
		System.out.println("  --ignore_type               - allow variants with different types to be merged");
		System.out.println("  --dup_to_ins                - convert duplications to insertions for SV merging and then convert them back");
		System.out.println("  --mark_specific             - mark calls in the original VCF files that have enough support to called specific");
		System.out.println("  --run_iris                  - run Iris before merging for refining the sequences of insertions");
		System.out.println();
		System.out.println("Notes:");
		System.out.println("  genome_file is required if the dup_to_ins option or the run_iris option is used.");
		System.out.println("  bam_list is required if the run_iris option is used.");
		System.out.println();
		
	}
	
	/*
	 * Method for parsing long integers allowing for abbreviations like 10k, 5m, etc.
	 */
	static long parseLong(String s) throws Exception
	{
		s = s.toLowerCase();
		if(s.endsWith("g") || s.endsWith("b") || s.endsWith("kkk"))
		{
			return (long)(Double.parseDouble(s.substring(0, s.length()-1)) * 1e9 + .5);
		}
		if(s.endsWith("m") || s.endsWith("kk"))
		{
			return (long)(Double.parseDouble(s.substring(0, s.length()-1)) * 1e6 + .5); 
		}
		if(s.endsWith("k"))
		{
			return (long)(Double.parseDouble(s.substring(0, s.length()-1)) * 1e3 + .5); 
		}
		return Long.parseLong(s);
	}
	
	static int parseInt(String s) throws Exception
	{
		return (int)parseLong(s);
	}
	
	/*
	 * Method for parsing all command line arguments and making sure the required
	 * arguments are given values
	 */
	static void parseArgs(String[] args) throws Exception
	{
		for(int i = 0; i<args.length; i++)
		{
			if(args[i].indexOf('=') == -1)
			{
				if(args[i].endsWith("ignore_strand"))
				{
					USE_STRAND = false;
				}
				else if(args[i].endsWith("ignore_type"))
				{
					USE_TYPE = false;
				}
				else if(args[i].endsWith("dup_to_ins"))
				{
					CONVERT_DUPLICATIONS = true;
				}
				else if(args[i].endsWith("mark_specific"))
				{
					MARK_SPECIFIC = true;
				}
				else if(args[i].endsWith("run_iris"))
				{
					RUN_IRIS = true;
				}
				continue;
			}
			int equalIdx = args[i].indexOf('=');
			String key = args[i].substring(0, equalIdx);
			while(key.length() > 0 && key.charAt(0) == '-')
			{
				key = key.substring(1);
			}
			String val = args[i].substring(1 + equalIdx);
			switch(key) 
			{
				case "max_dist":
					MAX_DIST = parseInt(val);
					break;
				case "max_dist_linear":
					MAX_DIST_LINEAR = Double.parseDouble(val);
					break;
				case "min_seq_id":
					MIN_SEQUENCE_SIMILARITY = Double.parseDouble(val);
					break;
				case "min_support":
					MIN_SUPPORT = parseInt(val);
					break;
				case "threads":
					THREADS = parseInt(val);
					break;
				case "spec_reads":
					SPECIFIC_MIN_RCOUNT = parseInt(val);
					break;
				case "spec_len":
					SPECIFIC_MIN_LENGTH = parseInt(val);
					break;
				case "file_list":
					FILE_LIST = val;
					break;
				case "out_file":
					OUT_FILE = val;
					break;
				case "genome_file":
					GENOME_FILE = val;
					break;
				case "bam_list":
					BAM_FILE_LIST = val;
					break;
				case "iris_args":
					IRIS_ARGS = val;
					break;
				default:
					break;
			}
		}
		if(FILE_LIST.length() == 0 || OUT_FILE.length() == 0)
		{
			usage();
			System.exit(1);
		}
		if(GENOME_FILE.length() == 0 && (RUN_IRIS || CONVERT_DUPLICATIONS))
		{
			usage();
			System.exit(1);
		}
		if(BAM_FILE_LIST.length() == 0 && RUN_IRIS)
		{
			usage();
			System.exit(1);
		}
		
		Path currentRelativePath = Paths.get("");
		OUT_DIR = currentRelativePath.toAbsolutePath().toString() + "/" + OUT_DIR;
		File f = new File(OUT_DIR);
		f.mkdir();
		
	}
}
