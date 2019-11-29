/*
 * Main interface for Jasmine
 */
import java.util.ArrayList;
import java.util.TreeMap;

public class Main {
public static void main(String[] args) throws Exception
{
	Settings.parseArgs(args);
	
	// The input file to SV merging may change based on the steps the user wants to run
	String currentInputFile = Settings.FILE_LIST;
	
	// Convert the duplications to insertions if the user wants to
	if(Settings.CONVERT_DUPLICATIONS)
	{
		currentInputFile = PipelineManager.convertDuplicationsToInsertions(currentInputFile);
	}
	
	// Mark calls with strong read support and long length as specific calls
	if(Settings.MARK_SPECIFIC)
	{
		currentInputFile = PipelineManager.markSpecificCalls(currentInputFile);
	}
	
	// Run iris if the user specifies that they want to run it
	if(Settings.RUN_IRIS)
	{
		currentInputFile = PipelineManager.runIris(currentInputFile);
	}
	
	// If running only preprocessing, stop now
	if(Settings.PREPROCESS_ONLY)
	{
		return;
	}
	
	// Get the variants and bin them into individual graphs
	TreeMap<String, ArrayList<Variant>> allVariants = VariantInput.readAllFiles(currentInputFile);
	
	// Initialize data structure for outputting merged variants
	VariantOutput output = new VariantOutput();
	
	// Get the number of samples to know the length of the SUPP_VEC field
	int sampleCount = VariantInput.countFiles(currentInputFile);
	
	// Merge each graph in parallel
	ParallelMerger pm = new ParallelMerger(allVariants, output, sampleCount);
	pm.run();
	
	System.out.println("Merging complete - outputting results");
	
	// Print the merged variants to a file if they have enough support
	output.writeMergedVariants(currentInputFile, Settings.OUT_FILE);
	
	System.out.println("Number of sets with multiple variants: " + pm.totalMerged.get()); 
	
	// Convert insertions back to duplications as needed
	if(Settings.CONVERT_DUPLICATIONS)
	{
		PipelineManager.convertInsertionsBackToDuplications();
	}
}
}
