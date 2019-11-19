/*
 * Utility for handling command line parameters and the usage message 
 */

public class Settings {
	
	static boolean USE_STRAND = true;
	static boolean USE_TYPE = true;
	static String FILE_LIST = "";
	static String OUT_FILE = "";
	static int MAX_DIST = 1000;
	static int MIN_SUPPORT = 2;
	
	static void usage()
	{
		System.out.println();
		System.out.println("Usage: java Main [args]");
		System.out.println("  Example: java Main file_list=filelist.txt out_file=out.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  file_list (String) - a file listing paths to all variant files to use (on separate lines)");
		System.out.println("  out_file  (String) - the name of the file to output the merged variants to");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  max_dist    (int) [1000] - the maximum distance variants can be apart when being merged");
		System.out.println("  min_support (int) [2]    - the maximum distance variants can be apart when being merged");
		System.out.println("  --ignore_strand          - allow variants with different strands to be merged");
		System.out.println("  --ignore_type            - allow variants with different types to be merged");
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
				case "min_support":
					MIN_SUPPORT = parseInt(val);
					break;
				case "file_list":
					FILE_LIST = val;
					break;
				case "out_file":
					OUT_FILE = val;
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
	}
}
