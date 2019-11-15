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
	int[] map; //negative if root, more negative means bigger set; if nonnegative, then it indicates the parent
	long[] sampleMask; // For each root, a bitmask of which samples are present in its component
	public Forest(Variant[] data)
	{
		int n = data.length;
		map = new int[n];
		Arrays.fill(map, -1);
		sampleMask = new long[n];
		for(int i = 0; i<n; i++) sampleMask[i] |= (1L << data[i].sample);
		
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
	public boolean union(int a, int b)
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
		if(map[roota] < map[rootb])
		{
			map[roota] += map[rootb]; //add the sizes
			map[rootb] = roota; //connect the smaller to the bigger
			sampleMask[roota] |= sampleMask[rootb];
		}
		else
		{
			map[rootb] += map[roota];
			map[roota] = rootb;
			sampleMask[rootb] |= sampleMask[roota];
		}
		return true;
	}
	
	private boolean okayEdge(int rootA, int rootB)
	{
		return (sampleMask[rootA] & sampleMask[rootB]) == 0;
	}
}
