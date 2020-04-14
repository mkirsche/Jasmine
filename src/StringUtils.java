import java.util.HashMap;

/*
 * A collection of String functions
 */
public class StringUtils {
	
	/*
	 * The sequence identity of two strings based on their edit distance
	 */
	static double editDistanceSimilarity(String s, String t)
	{
		int n = s.length(), m = t.length();
		int[][] editDistance = new int[n+1][m+1];
		for(int i = 1; i<=m; i++) editDistance[0][i] = i;
		for(int i = 1; i<=n; i++) editDistance[i][0] = i;
		for(int i = 1; i<=n; i++)
		{
			for(int j = 1; j<=m; j++)
			{
				boolean sameChar = s.charAt(i-1) == t.charAt(j-1);
				int bestDistance = editDistance[i-1][j-1] + (sameChar ? 0 : 1);
				bestDistance = Math.min(bestDistance, 1 + editDistance[i-1][j]);
				bestDistance = Math.min(bestDistance, 1 + editDistance[i][j-1]);
				editDistance[i][j] = bestDistance;
			}
		}
		return 1.0 * editDistance[n][m] / Math.max(n, m);
	}

	/*
	 * Gets the frequency of each k-mer in a string, skipping over non-base characters
	 */
	static HashMap<Integer, Integer> countKmers(String s, int k)
	{
		HashMap<Integer, Integer> kmerCount = new HashMap<Integer, Integer>();
		
		// The number of basepair characters (ACGT) we have seen so far
		int baseCount = 0;
		
		// The encoded (2 bits per character) value of the current kmer so far
		int kmer = 0;
		
		// Use a sliding window to get all of the kmer codes and add them to the frequency map
		for(int i = 0; i<s.length(); i++)
		{
			char c = s.charAt(i);
			
			// The value of the current character
			int base = -1;
			if(c == 'a' || c == 'A') base = 0;
			if(c == 'c' || c == 'C') base = 1;
			if(c == 'g' || c == 'G') base = 2;
			if(c == 't' || c == 'T') base = 3;
			
			// ONly consider basepair characters
			if(base != -1)
			{
				// Remove the base from the window that just left the sliding window and add the new base
				int allButTwoHighest = ((1 << (2*k-2)) - 1) & kmer;
				kmer = (allButTwoHighest << 2) | base;
				baseCount++;
				
				// If we have seen at least k bases, add the current kmer value to the map
				if(baseCount >= k)
				{
					if(kmerCount.containsKey(kmer))
					{
						kmerCount.put(kmer, kmerCount.get(kmer)+1);
					}
					else
					{
						kmerCount.put(kmer, 1);
					}
				}
			}
		}
		return kmerCount;
	}
	
	/*
	 * The sequence identity of two strings based on their kmer Jaccard distance
	 */
	static double jaccardSimilarity(String s, String t)
	{
		int k = Settings.K_JACCARD;
		
		// Get the frequencies of kmers in s
		HashMap<Integer, Integer> sKmerFreq = countKmers(s, k);
		if(sKmerFreq.size() <= 0)
		{
			return 1.0;
		}
		
		// Get the frequencies of kmer in t
		HashMap<Integer, Integer> tKmerFreq = countKmers(t, k);
		if(tKmerFreq.size() <= 0)
		{
			return 1.0;
		}
		
		// Compute the min and max count of each kmer to get the intersection and union, respectively 
		int intersection = 0, union = 0;
		
		// Iterate over everything in s - this includes both kmers distinct to s and kmers in both
		for(int sKmer : sKmerFreq.keySet())
		{
			int sFrequency = sKmerFreq.get(sKmer);
			int tFrequency = tKmerFreq.getOrDefault(sKmer, 0);
			intersection += Math.min(sFrequency,  tFrequency);
			union += Math.max(sFrequency,  tFrequency);
		}
		
		// Add the kmers unique to t to the union
		for(int tKmer : tKmerFreq.keySet())
		{
			if(!sKmerFreq.containsKey(tKmer))
			{
				union += tKmerFreq.get(tKmer);
			}
		}
		
		// Compute the Jaccard similarity as the intersection size divided by the union size
		return 1.0 * intersection / union;
		
	}
	
	/*
	 * Assumes input is a filename, and adds "_<desc>" right before the file extension
	 */
	static String addDescriptor(String input, String desc)
	{
		int idx = input.lastIndexOf(".");
		if(idx == -1)
		{
			return input + "_" + desc;
		}
		
		String before = input.substring(0, idx);
		String after = input.substring(idx);
		return before + "_" + desc + after;
	}
	
	/*
	 * Gets the basename of a file from its path by removing the directory name
	 */
	static String fileBaseName(String path)
	{
		int idx = path.lastIndexOf('/');
		if(idx == -1)
		{
			return path;
		}
		return path.substring(1 + idx);
	}
}
