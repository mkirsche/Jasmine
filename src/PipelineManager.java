import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Scanner;

public class PipelineManager {
	
/*
 * Convert duplications to insertions for all VCF files and update filelist
 * Returns a path to the new filelist
 */
static String convertDuplicationsToInsertions(String fileList) throws Exception
{
	String newFileList = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(fileList), "dupToIns");
	
	Scanner vcfListInput = new Scanner(new FileInputStream(new File(fileList)));
	ArrayList<String> vcfFiles = new ArrayList<String>();
			
	// Get a list of all VCF Files to refine
	while(vcfListInput.hasNext())
	{
		String line = vcfListInput.nextLine();
		if(line.length() > 0)
		{
			vcfFiles.add(line);
		}
	}
	
	PrintWriter newFileListOut = new PrintWriter(new File(newFileList));
	
	for(int i = 0; i<vcfFiles.size(); i++)
	{
		String vcfFile = vcfFiles.get(i);
		String newVcfFile = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(vcfFile), "dupToIns");
		newFileListOut.println(newVcfFile);
		DuplicationsToInsertions.convertFile(vcfFile, Settings.GENOME_FILE, newVcfFile);
	}
	vcfListInput.close();
	newFileListOut.close();
		
	return newFileList;
}

static String runIris(String fileList) throws Exception
{
	String refinedInput = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(fileList), "irisRefined");
	Scanner vcfListInput = new Scanner(new FileInputStream(new File(fileList)));
	Scanner bamListInput = new Scanner(new FileInputStream(new File(Settings.BAM_FILE_LIST)));
	ArrayList<String> vcfFiles = new ArrayList<String>(), bamFiles = new ArrayList<String>();
	
	PrintWriter newFileListOut = new PrintWriter(new File(refinedInput));
	
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
		String newVcfFile = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(vcfFile), "irisRefined");
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
	return refinedInput;
}

/*
 * Converts insertions in the output file which used to include a duplication back to their original types
 * Moves the old output file and replaces it with the updated one
 */
static void convertInsertionsBackToDuplications() throws Exception
{
	String unconvertedOutput = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(Settings.OUT_FILE), "dupToIns");
	File f;
	if((f = new File(unconvertedOutput)).exists())
	{
		f.delete();
	}
	Files.move(Paths.get(Settings.OUT_FILE), Paths.get(unconvertedOutput));
	InsertionsToDuplications.convertFile(unconvertedOutput, Settings.OUT_FILE);
}

}
