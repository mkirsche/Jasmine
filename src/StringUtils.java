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
	 * The sequence identity of two strings based on their kmer Jaccard distance
	 */
	static int k = 9;
	static double jaccardSimilarity(String s, String t)
	{
		HashMap<Integer, Integer> kmerCount = new HashMap<Integer, Integer>();
		int baseCount = 0;
		int n = s.length(), m = t.length();
		int kmer = 0;
		for(int i = 0; i<n; i++)
		{
			char c = s.charAt(i);
			int base = -1;
			if(c == 'a' || c == 'A') base = 0;
			if(c == 'c' || c == 'C') base = 1;
			if(c == 'g' || c == 'G') base = 2;
			if(c == 't' || c == 'T') base = 3;
			if(base != -1)
			{
				int allButTwoHighest = ((1 << (2*k-2)) - 1) & kmer;
				kmer = (allButTwoHighest << 2) | base;
				baseCount++;
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
		
		int totalS = baseCount - (k - 1);
		if(totalS <= 0)
		{
			return 1.0;
		}
		
		baseCount = 0;
		int intersect = 0;
		for(int i = 0; i<m; i++)
		{
			char c = t.charAt(i);
			int base = -1;
			if(c == 'a' || c == 'A') base = 0;
			if(c == 'c' || c == 'C') base = 1;
			if(c == 'g' || c == 'G') base = 2;
			if(c == 't' || c == 'T') base = 3;
			if(base != -1)
			{
				int allButTwoHighest = ((1 << (2*k-2)) - 1) & kmer;
				kmer = (allButTwoHighest << 2) | base;
				baseCount++;
				if(baseCount >= k)
				{
					if(kmerCount.containsKey(kmer))
					{
						intersect++;
						kmerCount.put(kmer, kmerCount.get(kmer)-1);
						if(kmerCount.get(kmer) == 0)
						{
							kmerCount.remove(kmer);
						}
					}
				}
			}
		}
		
		int totalT = baseCount - (k - 1);
		
		if(totalT <= 0)
		{
			return 1.0;
		}
		return 1.0 * intersect / Math.max(totalS, totalT);
		
	}
	
	/*
	 * Assumes input is a filename, and adds "_<desc>: right before the file extension
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
	 * Gets the name of a file from its path by removing the directory name
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
