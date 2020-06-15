/*
 * A utility for managing pipeline steps for multiple VCF files
 * Most of the pre-processing and post-processing steps are done on a per-VCF basis,
 * so this manager performs them for all files and updates the filelist to point to a
 * list of updated files instead of the original ones
 */

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
	ArrayList<String> vcfFiles = getFilesFromList(fileList);
	ArrayList<String> newVcfFiles = new ArrayList<String>();	
	
	for(int i = 0; i<vcfFiles.size(); i++)
	{
		String vcfFile = vcfFiles.get(i);
		String newVcfFile = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(vcfFile), "dupToIns");
		newVcfFiles.add(newVcfFile);
		DuplicationsToInsertions.convertFile(vcfFile, Settings.GENOME_FILE, newVcfFile);
	}
		
	return buildUpdatedFileList(fileList, "dupToIns", newVcfFiles);
}

/*
 * Normalize types for all VCF files and update filelist
 * Returns a path to the new filelist
 */
static String normalizeTypes(String fileList) throws Exception
{
	
	ArrayList<String> vcfFiles = getFilesFromList(fileList);
	ArrayList<String> newVcfFiles = new ArrayList<String>();
	
	for(int i = 0; i<vcfFiles.size(); i++)
	{
		String vcfFile = vcfFiles.get(i);
		
		String newVcfFile = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(vcfFile), "normalizeTypes");
		newVcfFiles.add(newVcfFile);
		
		NormalizeTypes.convertFile(vcfFile, newVcfFile);
	}
	
	return buildUpdatedFileList(fileList, "normalizeTypes", newVcfFiles);
	
}

/*
 * Run Iris on all VCF files and update the filelist
 * Returns a path to the new filelist
 */
static String runIris(String fileList) throws Exception
{
	ArrayList<String> vcfFiles = getFilesFromList(fileList), bamFiles = getFilesFromList(Settings.BAM_FILE_LIST);
	ArrayList<String> newVcfFiles = new ArrayList<String>();
	
	// Get any optional arguments to be passed to Iris that the user specified
	String[] optionalArgs = Settings.IRIS_ARGS.split(",");
	
	// Refine one VCF file at a time
	for(int i = 0; i<vcfFiles.size(); i++)
	{
		String vcfFile = vcfFiles.get(i);
		String newVcfFile = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(vcfFile), "irisRefined");
		String bamFile = bamFiles.get(i);
		newVcfFiles.add(newVcfFile);
		
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
	
	return buildUpdatedFileList(fileList, "irisRefined", newVcfFiles);
}

/*
 * Mark all specific calls in the input VCFs and update the filelist
 * Returns a path to the new filelist
 */
static String markSpecificCalls(String fileList) throws Exception
{
	ArrayList<String> vcfFiles = getFilesFromList(fileList);
	ArrayList<String> newVcfFiles = new ArrayList<String>();
			
	//PrintWriter newFileListOut = new PrintWriter(new File(newFileList));
	
	for(int i = 0; i<vcfFiles.size(); i++)
	{
		String vcfFile = vcfFiles.get(i);
		String newVcfFile = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(vcfFile), "markedSpec");
		newVcfFiles.add(newVcfFile);
		MarkSpecificCalls.convertFile(vcfFile, newVcfFile, Settings.SPECIFIC_MIN_RCOUNT, Settings.SPECIFIC_MIN_LENGTH);
	}
	
	return buildUpdatedFileList(fileList, "markedSpec", newVcfFiles);
}

/*
 * Builds an updated file list after running a pre-processing step
 */
static String buildUpdatedFileList(String oldFileList, String suffix, ArrayList<String> newVcfFiles) throws Exception
{
	if(Settings.USING_FILE_LIST)
	{
		String newFileList = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(oldFileList), suffix);
		PrintWriter newFileListOut = new PrintWriter(new File(newFileList));
		for(String newVcfFile : newVcfFiles)
		{
			newFileListOut.println(newVcfFile);
		}
		newFileListOut.close();
		return newFileList;
	}
	else
	{
		StringBuilder res = new StringBuilder("");
		for(int i = 0; i<newVcfFiles.size(); i++)
		{
			res.append(newVcfFiles.get(i));
			if(i < newVcfFiles.size() - 1)
			{
				res.append(",");
			}
		}
		return res.toString();
	}
}

/*
 * Reverts insertions in the output which were formerly duplications back to duplication calls 
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

/*
 * Adds genotypes to the merged VCF based on the genotypes in individual samples
 * Moves the old output file and replaces it with the updated one
 */
static void addGenotypes(String fileList) throws Exception
{
	String unconvertedOutput = Settings.OUT_DIR + "/" + StringUtils.addDescriptor(StringUtils.fileBaseName(Settings.OUT_FILE), "noGenotypes");
	File f;
	if((f = new File(unconvertedOutput)).exists())
	{
		f.delete();
	}
	Files.move(Paths.get(Settings.OUT_FILE), Paths.get(unconvertedOutput));
	AddGenotypes.addGenotypes(unconvertedOutput, fileList, Settings.OUT_FILE);
}

/*
 * Reads the list of files from either a specified list file or the comma-separated command line argument
 */
static ArrayList<String> getFilesFromList(String fileList) throws Exception
{
	ArrayList<String> res = new ArrayList<String>();
	
	if(!Settings.USING_FILE_LIST)
	{
		String[] fns = fileList.split(",");
		for(String fn : fns) res.add(fn);
		return res;
	}
	
	if(new File(fileList).exists())
	{
		Scanner vcfListInput = new Scanner(new FileInputStream(new File(fileList)));
				
		while(vcfListInput.hasNext())
		{
			String line = vcfListInput.nextLine();
			if(line.length() > 0)
			{
				res.add(line);
			}
		}
		vcfListInput.close();
	}
	
	return res;
}

}
