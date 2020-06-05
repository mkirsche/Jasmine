/*
 * A representation of a forest using a union-find data structure
 * It allows nodes to be merged, checking if their components share
 * any variants from the same sample.  For now, only up to 64 samples
 * are supported, but sampleMask can be replaced with actual bitsets
 * at a small cost to runtime.
 */

import java.util.Arrays;

public class Forest
{
	int[] map; // map[i] is negative if root, more negative means bigger set; if nonnegative, then it indicates the parent
	long[][] sampleMask; // For each root node, a bitmask of which samples are present in its component
	static int samplesPerMask = 63;
	
	public Forest(Variant[] data)
	{
		int n = data.length;
		int maxSample = 0;
		for(int i = 0; i<n; i++)
		{
			maxSample = Math.max(maxSample, data[i].sample);
		}
		
		// Each component may require multiple 64-bit integers to hold its bitset of sample IDs if there are many samples
		int masksNeeded = maxSample / samplesPerMask + 1;
		map = new int[n];
		Arrays.fill(map, -1);
		sampleMask = new long[masksNeeded][n];
		for(int i = 0; i<n; i++)
		{
			int maskId = data[i].sample / samplesPerMask;
			int maskVal = data[i].sample % samplesPerMask;
			sampleMask[maskId][i] |= (1L << maskVal);
		}
		
	}
	
	/*
	 * Get the root of the component containing a variant
	 */
	public int find(int x)
	{
		if(map[x] < 0)
			return x;
		else
		{
			map[x] = find(map[x]);
			return map[x];
		}
	}
	
	/*
	 * Add an edge between two variants
	 */
	public boolean canUnion(int a, int b)
	{
		int roota = find(a), rootb = find(b);
		if(roota == rootb)
		{
			return false;
		}
		if(!okayEdge(roota, rootb))
		{
			return false;
		}
		
		return true;
	}
	
	public void union(int a, int b)
	{
		int roota = find(a), rootb = find(b);
		if(map[roota] < map[rootb])
		{
			map[roota] += map[rootb]; //add the sizes
			map[rootb] = roota; //connect the smaller to the bigger
			for(int j = 0; j<sampleMask.length; j++)
			{
				sampleMask[j][roota] |= sampleMask[j][rootb];
			}
		}
		else
		{
			map[rootb] += map[roota];
			map[roota] = rootb;
			for(int j = 0; j<sampleMask.length; j++)
			{
				sampleMask[j][rootb] |= sampleMask[j][roota];
			}
		}
	}
	
	/*
	 * Whether or not adding an edge will avoid causing intrasample merging
	 */
	private boolean okayEdge(int rootA, int rootB)
	{
		if(Settings.ALLOW_INTRASAMPLE)
		{
			return true;
		}
		for(int j = 0; j<sampleMask.length; j++)
		{
			if((sampleMask[j][rootA] & sampleMask[j][rootB]) != 0)
			{
				return false;
			}
		}
		return true;
	}
}
