/*
 * Main interface for Thriver.
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.TreeMap;

public class Main {
public static void main(String[] args) throws Exception
{
	Settings.parseArgs(args);
	
	/*
	 * Run pre-processing steps as needed
	 */
	
	// The input file to SV merging may change based on the steps the user wants to run
	String currentInputFile = Settings.FILE_LIST;
	
	/*
	 * Convert the duplications to insertions if the user wants to
	 */
	if(Settings.CONVERT_DUPLICATIONS)
	{
		String typeConvertedInput = StringUtils.addDescriptor(currentInputFile, "dupToIns");
		DuplicationsToInsertions.convertFile(currentInputFile, Settings.GENOME_FILE, typeConvertedInput);
		currentInputFile = typeConvertedInput; 
	}
	
	/*
	 * Run iris if the user specifies they want it run
	 */
	if(Settings.RUN_IRIS)
	{
		String refinedInput = StringUtils.addDescriptor(currentInputFile, "irisRefined");
		Scanner vcfListInput = new Scanner(new FileInputStream(new File(currentInputFile)));
		Scanner bamListInput = new Scanner(new FileInputStream(new File(Settings.BAM_FILE_LIST)));
		ArrayList<String> vcfFiles = new ArrayList<String>(), bamFiles = new ArrayList<String>();
		
		PrintWriter newFileListOut =new PrintWriter(new File(refinedInput));
		
		// Get a list of all VCF Files to refine
		while(vcfListInput.hasNext())
		{
			String line = vcfListInput.nextLine();
			if(line.length() > 0)
			{
				vcfFiles.add(line);
			}
		}
		
		// Get a list of the corresponding read alignment files for each VCF to refine
		while(bamListInput.hasNext())
		{
			String line = bamListInput.nextLine();
			if(line.length() > 0)
			{
				bamFiles.add(line);
			}
		}
		
		// Get any optional arguments to be passed to Iris that the user specified
		String[] optionalArgs = Settings.IRIS_ARGS.split(",");
		
		// Refine one VCF file at a time
		for(int i = 0; i<vcfFiles.size(); i++)
		{
			String vcfFile = vcfFiles.get(i);
			String newVcfFile = StringUtils.addDescriptor(vcfFile, "irisRefined");
			String bamFile = bamFiles.get(i);
			newFileListOut.println(newVcfFile);
			
			// Generate the required Iris arguments based on the filenames
			String[] requiredArgs = new String[]
			{
				"genome_in=" + Settings.GENOME_FILE,
				"vcf_in=" + vcfFile,
				"reads_in=" + bamFile,
				"vcf_out=" + newVcfFile
			};
			
			// Concatenate the required arguments and the optional ones
			String[] allArgs = new String[optionalArgs.length + requiredArgs.length];
			for(int j = 0; j<requiredArgs.length; j++)
			{
				allArgs[j] = requiredArgs[j];
			}
			for(int j = 0; j<optionalArgs.length; j++)
			{
				allArgs[j + requiredArgs.length] = optionalArgs[j];
			}
			
			// Actually run Iris
			Iris.runIris(allArgs);
		}
		
		// Close all Scanners and PrintWriters
		vcfListInput.close();
		bamListInput.close();
		newFileListOut.close();
		
		// Update the input filename to be the refined one
		currentInputFile = refinedInput;
	}
	
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
	
	/*
	 * Run post-processing steps as needed
	 */
	if(Settings.CONVERT_DUPLICATIONS)
	{
		String unconvertedOutput = StringUtils.addDescriptor(Settings.OUT_FILE, "dupToIns");
		Files.move(Paths.get(Settings.OUT_FILE), Paths.get(unconvertedOutput));
		InsertionsToDuplications.convertFile(unconvertedOutput, Settings.OUT_FILE);
	}
}
}
