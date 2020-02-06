/*
 * Utility for handling command line parameters and the usage message 
 */

import java.io.File;
import java.nio.file.Path;
import java.nio.file.Paths;

public class Settings {

	/*
	 * Command line parameters - see the usage menu for what can be changed and how
	 */
	static boolean USE_STRAND = true;
	static boolean USE_TYPE = true;
	static String FILE_LIST = "";
	static String OUT_FILE = "";
	static int MAX_DIST = 1000;
	static double MAX_DIST_LINEAR = 0.0;
	static int MIN_SUPPORT = 1;
	static double MIN_SEQUENCE_SIMILARITY = 0;
	static boolean USE_EDIT_DISTANCE = false;
	static int K_JACCARD = 9;
	static int MAX_DUP_LEN = 10000;
	static int KD_TREE_NORM = 2;
	static boolean CHANGE_VAR_IDS = true;
	static boolean USE_END = false;
	static boolean MAX_DIST_SET = false;
	static int MIN_DIST = -1; // -1 means no minimum
	static boolean OUTPUT_GENOTYPES = false;
	static boolean INPUTS_MERGED = true;
	
	static String SAMTOOLS_PATH = "samtools";
	
	static boolean PREPROCESS_ONLY = false;
	static boolean POSTPROCESS_ONLY = false;
	static boolean CONVERT_DUPLICATIONS = false;
	static boolean MARK_SPECIFIC = false;
	static boolean RUN_IRIS = false;
	static boolean FIX_ENDS = true;
	static String GENOME_FILE = "";
	static String BAM_FILE_LIST = "";
	static String IRIS_ARGS = "";
	
	static String OUT_DIR = "output";
	static int THREADS = 1;
	
	static int SPECIFIC_MIN_RCOUNT = 10;
	static int SPECIFIC_MIN_LENGTH = 30;
	
	static boolean CENTROID_MERGE = false;
	static boolean CLIQUE_MERGE = false;
	
	static boolean ALLOW_INTRASAMPLE = false;
	static boolean NORMALIZE_TYPE = false;
	
	/*
	 * Print the usage menu
	 */
	static void usage()
	{
		System.out.println();
		System.out.println("Jasmine version 1.0.1");
		System.out.println("Usage: java -cp src Main [args]");
		System.out.println("  Example: java -cp src Main file_list=filelist.txt out_file=out.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  file_list (String) - a file listing paths to all variant files to use (on separate lines)");
		System.out.println("  out_file  (String) - the name of the file to output the merged variants to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  max_dist        (int)    [1000]     - the maximum distance variants can be apart when being merged");
		System.out.println("  min_dist        (int)    [-1]       - the minimum distance threshold a variant can have when using max_dist_linear");
		System.out.println("  max_dist_linear (float)  [0]        - make max_dist this proportion of the length of each variant (overrides max_dost)");
		System.out.println("  kd_tree_norm    (int)    [2]        - the power to use in kd-tree distances (1 is Manhattan, 2 is Euclidean, etc.)");
		System.out.println("  min_seq_id      (float)  [0]        - the minimum sequence identity for two insertions to be merged");
		System.out.println("  k_jaccard       (int)    [9]        - the kmer size to use when computing Jaccard similarity of insertions");
		System.out.println("  max_dup_length  (int)    [10k]      - the maximum length of duplication that can be converted to an insertion");
		System.out.println("  min_support     (int)    [1]        - the minimum number of callsets a variant must be in to be output");
		System.out.println("  threads         (int)    [2]        - the number of threads to use for merging the variants");
		System.out.println("  spec_reads      (int)    [10]       - the minimum number of reads a variant needs to be in the specific callset");
		System.out.println("  spec_len        (int)    [30]       - the minimum length a variant needs to be in the specific callset");
		System.out.println("  genome_file     (String) []         - the reference genome being used");
		System.out.println("  bam_list        (String) []         - a file listing paths to BAMs in the same order as the VCFs");
		System.out.println("  iris_args       (String) []         - a comma-separated list of optional arguments to pass to Iris");
		System.out.println("  out_dir         (String) [output]   - the directory where intermediate files go");
		System.out.println("  samtools_path   (String) [samtools] - the path to the samtools executable used for coverting duplications");
		System.out.println("  --ignore_strand                     - allow variants with different strands to be merged");
		System.out.println("  --ignore_type                       - allow variants with different types to be merged");
		System.out.println("  --dup_to_ins                        - convert duplications to insertions for SV merging and then convert them back");
		System.out.println("  --mark_specific                     - mark calls in the original VCF files that have enough support to called specific");
		System.out.println("  --run_iris                          - run Iris before merging for refining the sequences of insertions");
		System.out.println("  --use_edit_dist                     - use edit distance for comparing insertion sequences instead of Jaccard");
		System.out.println("  --preprocess_only                   - only run the preprocessing and not the actual merging or post-processing");
		System.out.println("  --postprocess_only                  - only run the postprocessing and not the actual merging or pre-processing");
		System.out.println("  --keep_var_ids                      - don't change variant IDs (should only be used if input IDs are unique across samples)");
		System.out.println("  --use_end                           - use the end coordinate as the second coordinate instead of the variant length");
		System.out.println("  --output_genotypes                  - print the genotypes of the consensus variants in all of the samples they came from");
		System.out.println("  --ignore_merged_inputs              - ignore merging info such as support vectors which is already present in the inputs");
		System.out.println("  --centroid_merging                  - require every group to have a centroid which is within the distance threshold of each variant");
		System.out.println("  --clique_merging                    - require every group to have each pair within in it be mergeable");
		System.out.println("  --allow_intrasample                 - allow variants in the same sample to be merged");
		System.out.println("  --normalize_type                    - convert all variants to INS/DEL/DUP/INV/TRA");
		System.out.println("  --leave_breakpoints                 - leave breakpoints as they are even if they are inconsistent");
		System.out.println();
		System.out.println("Notes:");
		System.out.println("  genome_file is required if the dup_to_ins option or the run_iris option is used.");
		System.out.println("  bam_list is required if the run_iris option is used.");
		System.out.println("  Setting both max_dist_linear and max_dist sets thresholds to minimum of max_dist and max_dist_linear * sv_length");
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
		if(args.length == 1 && (args[0].equalsIgnoreCase("--version") || args[0].equalsIgnoreCase("-v")))
		{
			System.out.println("Jasmine version 1.0.0");
			System.exit(0);
		}
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
				else if(args[i].endsWith("use_edit_dist"))
				{
					USE_EDIT_DISTANCE = true;
				}
				else if(args[i].endsWith("preprocess_only"))
				{
					PREPROCESS_ONLY = true;
				}
				else if(args[i].endsWith("postprocess_only"))
				{
					POSTPROCESS_ONLY = true;
				}
				else if(args[i].endsWith("keep_var_ids"))
				{
					CHANGE_VAR_IDS = false;
				}
				else if(args[i].endsWith("use_end"))
				{
					USE_END = true;
				}
				else if(args[i].endsWith("output_genotypes"))
				{
					OUTPUT_GENOTYPES = true;
				}
				else if(args[i].endsWith("ignore_merged_inputs"))
				{
					INPUTS_MERGED = false;
				}
				else if(args[i].endsWith("centroid_merging"))
				{
					CENTROID_MERGE = true;
				}
				else if(args[i].endsWith("clique_merging"))
				{
					CLIQUE_MERGE = true;
				}
				else if(args[i].endsWith("allow_intrasample"))
				{
					ALLOW_INTRASAMPLE = true;
				}
				else if(args[i].endsWith("normalize_type"))
				{
					NORMALIZE_TYPE = true;
				}
				else if(args[i].endsWith("leave_breakpoints"))
				{
					FIX_ENDS = false;
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
					MAX_DIST_SET = true;
					break;
				case "min_dist":
					MIN_DIST = parseInt(val);
					break;
				case "max_dist_linear":
					MAX_DIST_LINEAR = Double.parseDouble(val);
					break;
				case "kd_tree_norm":
					KD_TREE_NORM = parseInt(val);
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
				case "k_jaccard":
					K_JACCARD = parseInt(val);
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
				case "out_dir":
					OUT_DIR = val;
					break;
				case "samtools_path":
					SAMTOOLS_PATH = val;
					break;
				case "max_dup_length":
					MAX_DUP_LEN = parseInt(val);
					break;
				default:
					break;
			}
		}
		if(FILE_LIST.length() == 0 && !POSTPROCESS_ONLY)
		{
			usage();
			System.exit(1);
		}
		if(OUT_FILE.length() == 0 && !PREPROCESS_ONLY)
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
