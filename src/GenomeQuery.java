/*
 * A wrapper around samtools to get substrings of a genome
 */

import java.io.File;
import java.io.InputStream;
import java.util.Scanner;

public class GenomeQuery {
	String filename;
	
	/*
	 * The constructor tests that samtools works and that the genome file exists
	 */
	GenomeQuery(String filename) throws Exception
	{
		testSamtoolsInstalled();
		boolean validFile = new File(filename).exists();
		if(!validFile)
		{
			throw new Exception("geonome file does not exist: " + filename);
		}
		this.filename = filename;
	}
	
	/*
	 * Runs a simple samtools command and inspects the exit code to make sure it is installed
	 */
	void testSamtoolsInstalled() throws Exception
	{
		String samtoolsTestCommand = Settings.SAMTOOLS_PATH;
		Process child = Runtime.getRuntime().exec(samtoolsTestCommand);
        int seqExit = child.waitFor();
		
        // Exit code > 1 means the command failed, usually because samtools is not installed or on path
        if(seqExit > 1)
        {
        	throw new Exception("samtools produced bad exit code (" + seqExit + ") - check path: " + Settings.SAMTOOLS_PATH);
        }
	}
	
	/*
	 * Queries a genomic substring - runs samtools faidx <genomeFile> chr:startPos-endPos
	 */
	String genomeSubstring(String chr, long startPos, long endPos) throws Exception
	{
		if(startPos > endPos)
		{
			return "";
		}
		String faidxCommand = String.format("%s faidx %s %s:%d-%d", Settings.SAMTOOLS_PATH, filename, chr, startPos, endPos);
		Process child = Runtime.getRuntime().exec(faidxCommand);
        InputStream seqStream = child.getInputStream();
		Scanner seqInput = new Scanner(seqStream);
		
		// Make sure it produced an actual output
		if(!seqInput.hasNext())
        {
        	seqInput.close();
        	throw new Exception("samtools faidx did not produce an output: " + faidxCommand);
        }
		// Read in and ignore sequence name
		seqInput.next();
		
		// Make sure there's a sequence
		if(!seqInput.hasNext())
		{
			seqInput.close();
        	throw new Exception("samtools faidx produced a sequence name but not an actual sequence: " + faidxCommand);
		}
		
		// Concatenate all lines of the output sequence
		StringBuilder res = new StringBuilder("");
		while(seqInput.hasNext())
		{
			res.append(seqInput.next());
		}
		seqInput.close();

		return res.toString();
	}
}
