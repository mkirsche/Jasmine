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
}
