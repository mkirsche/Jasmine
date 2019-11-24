/*
 * Main interface for Thriver.
 */
import java.util.ArrayList;
import java.util.TreeMap;

public class Main {
public static void main(String[] args) throws Exception
{
	Settings.parseArgs(args);
	
	// The input file to SV merging may change based on the steps the user wants to run
	String currentInputFile = Settings.FILE_LIST;
	
	/*
	 * Convert the duplications to insertions if the user wants to
	 */
	if(Settings.CONVERT_DUPLICATIONS)
	{
		currentInputFile = PipelineManager.convertDuplicationsToInsertions(currentInputFile);
	}
	
	if(Settings.MARK_SPECIFIC)
	{
		currentInputFile = PipelineManager.markSpecificCalls(currentInputFile);
	}
	
	/*
	 * Run iris if the user specifies they want it run
	 */
	if(Settings.RUN_IRIS)
	{
		currentInputFile = PipelineManager.runIris(currentInputFile);
	}
	
	// Get the variants and bin them into individual graphs
	TreeMap<String, ArrayList<Variant>> allVariants = VariantInput.readAllFiles(currentInputFile);
	
	VariantOutput output = new VariantOutput();
	
	// Get the number of samples to know the length of the SUPP_VEC field
	int sampleCount = VariantInput.countFiles(currentInputFile);
	
	// Merge each graph in parallel
	ParallelMerger pm = new ParallelMerger(allVariants, output, sampleCount);
	pm.run();
	
	// Print the merged variants to a file if they have enough support
	output.writeMergedVariants(currentInputFile, Settings.OUT_FILE, Settings.MIN_SUPPORT);
	System.out.println("Number of sets with multiple variants: " + pm.totalMerged.get()); 
	
	// Convert insertions back to duplications as needed
	if(Settings.CONVERT_DUPLICATIONS)
	{
		PipelineManager.convertInsertionsBackToDuplications();
	}
}
}
