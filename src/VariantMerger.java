/*
 * Main interface for merging variants.  This assumes all variants are on the
 * same chromosome, and have the same type/strand if that separation is desired.
 */

import java.util.ArrayList;
import java.util.PriorityQueue;

public class VariantMerger
{
	// An array of all of the variants to be considered
	Variant[] data;
	
	// The number of total variants
	int n;
	
	// A forest in which connected components will represent merged groups
	Forest forest;
	
	// A KD-tree data structure for fast k-nearest-neighbors queries
	KDTree knn;
	
	// Indices of variants in each group, used for more advanced distance checks like clique and centroid
	ArrayList<Integer>[] merged;

	@SuppressWarnings("unchecked")
	public VariantMerger(Variant[] data)
	{
		n = data.length;
		
		forest = new Forest(data);
		knn = new KDTree(data);
		
		for(int i = 0; i<n; i++) data[i].index = i;
		this.data = data;
		
		if(Settings.CENTROID_MERGE || Settings.CLIQUE_MERGE)
		{
			merged = new ArrayList[n];
			for(int i = 0; i<n; i++)
			{
				merged[i] = new ArrayList<Integer>();
				merged[i].add(i);
			}
		}
	}
	
	/*
	 * Helper function to convert an ArrayList to an array to make the constructor more flexible
	 */
	static Variant[] listToArray(ArrayList<Variant> data)
	{
		int length = data.size();
		Variant[] asArray = new Variant[length];
		for(int i = 0; i<length; i++)
		{
			asArray[i] = data.get(i);
		}
		return asArray;
	}
	
	/*
	 * Alternate constructor which takes a list instead of an array
	 */
	public VariantMerger(ArrayList<Variant> data)
	{
		this(listToArray(data));
	}
	
	/*
	 * Runs the core algorithm for building the implicit merging graph and
	 * performing merging
	 */
	void runMerging()
	{
		if(n == 1)
		{
			return;
		}
		
		// For each variant v, how many of its nearest neighbors have had their edges
		// from v considered already.  
		int[] countEdgesProcessed = new int[n];
		
		// nearestNeighbors will be used as a cache to store the next few nearest neighbors
		// The purpose of this is to prevent performing a new KNN-query every time an edge
		// is considered, but instead a logarithmic number of times.
		Variant[][] nearestNeighbors = new Variant[n][];
		
		// A heap of edges to be processed in non-decreasing order of distance
		PriorityQueue<Edge> toProcess = new PriorityQueue<Edge>();
		
		// Get the first 4 nearest neighbors for every variant, and add their first edges to
		// the heap
		for(int i = 0; i<n; i++)
		{
			nearestNeighbors[i] = knn.kNearestNeighbor(data[i], 4);
			toProcess.add(new Edge(i, nearestNeighbors[i][0].index, data[i].distance(nearestNeighbors[i][0])));
			countEdgesProcessed[i]++;
		}
		
		while(!toProcess.isEmpty())
		{
			Edge e = toProcess.poll();
			boolean valid = forest.canUnion(e.from, e.to);
			if(valid)
			{
				int fromRoot = forest.find(e.from);
				int toRoot = forest.find(e.to);
				if(Settings.CLIQUE_MERGE || Settings.CENTROID_MERGE)
				{
					// Makes sure all newly merged variant pairs are within the distance threshold
					if(Settings.CLIQUE_MERGE)
					{
						for(int i = 0; i<merged[fromRoot].size() && valid; i++)
						{
							Variant candidateFrom = data[merged[fromRoot].get(i)];
							for(int j = 0; j<merged[toRoot].size() && valid; j++)
							{
								Variant candidateTo = data[merged[toRoot].get(j)];
								int maxDistAllowed = Math.max(candidateFrom.maxDist, candidateTo.maxDist);
								if(candidateFrom.distance(candidateTo) > maxDistAllowed + 1e-9)
								{
									valid = false;
								}
							}
						}
					}
					// Make sure everything being merged can be merged with their overall centroid
					else if(Settings.CENTROID_MERGE)
					{
						double avgStart = 0.0, avgEnd = 0.0;
						for(int i = 0; i<merged[fromRoot].size(); i++)
						{
							Variant v = data[merged[fromRoot].get(i)];
							avgStart += v.start;
							avgEnd += v.end;
						}
						for(int i = 0; i<merged[toRoot].size(); i++)
						{
							Variant v = data[merged[toRoot].get(i)];
							avgStart += v.start;
							avgEnd += v.end;
						}
						avgStart /= merged[fromRoot].size() + merged[toRoot].size();
						avgEnd /= merged[fromRoot].size() + merged[toRoot].size();
						
						for(int i = 0; i<merged[fromRoot].size() && valid; i++)
						{
							Variant v = data[merged[fromRoot].get(i)];
							valid &= v.distFromPoint(avgStart, avgEnd) <= v.maxDist + 1e-9;
						}
						for(int i = 0; i<merged[toRoot].size() && valid; i++)
						{
							Variant v = data[merged[toRoot].get(i)];
							valid &= v.distFromPoint(avgStart, avgEnd) <= v.maxDist + 1e-9;
						}
					}
					if(valid)
					{
						forest.union(fromRoot, toRoot);
						if(forest.map[fromRoot] < 0)
						{
							// The first variant is the new root of the union-find component
							for(int x : merged[toRoot])
							{
								merged[fromRoot].add(x);
							}
						}
						else
						{
							for(int x : merged[fromRoot])
							{
								merged[toRoot].add(x);
							}
						}
					}
				}
				
				// Two variants are being merged here - nothing needs to be done
				else
				{
					forest.union(e.from, e.to);
				}
			}
			
			while(true)
			{
				
				// If we already used the stored neighbors, query again for twice as many
				if(countEdgesProcessed[e.from] >= nearestNeighbors[e.from].length)
				{
					nearestNeighbors[e.from] = knn.kNearestNeighbor(data[e.from], 2 * nearestNeighbors[e.from].length);
				}
				
				// If we tried to get more and didn't find anymore, then we are done with this variant
				if(countEdgesProcessed[e.from] >= nearestNeighbors[e.from].length)
				{
					break;
				}
				Variant candidateTo = nearestNeighbors[e.from][countEdgesProcessed[e.from]];
				
				// This edge was invalid because of distance from the query, so stop looking at any edges 
				// since they'll only get farther away
				int maxDistAllowed = Math.max(data[e.from].maxDist, candidateTo.maxDist);
				if(data[e.from].distance(candidateTo) > maxDistAllowed + 1e-9)
				{
					break;
				}
				
				// If edge was invalid because of coming from the same sample, ignore it and try the next one
				else if(data[e.from].sample == candidateTo.sample)
				{
					countEdgesProcessed[e.from]++;
					continue;
				}
				
				// If sequences weren't similar enough for two insertions, ignore and try again
				else if(!data[e.from].passesStringSimilarity(candidateTo))
				{
					countEdgesProcessed[e.from]++;
					continue;
				}
				
				// The next edge is something we want to consider since it is close enough and goes to a
				// different sample
				else
				{
					toProcess.add(new Edge(e.from, candidateTo.index, data[e.from].distance(candidateTo)));
					countEdgesProcessed[e.from]++;
					break;
				}
			}
		}
	}
	
	/*
	 * Get an array of all of the groups of variants
	 */
	@SuppressWarnings("unchecked")
	ArrayList<Variant>[] getGroups()
	{
		ArrayList<Variant>[] res = new ArrayList[n];
		for(int i = 0; i<n; i++)
		{
			res[i] = new ArrayList<Variant>();
		}
		for(int i = 0; i<n; i++)
		{
			if(forest.map[i] < 0)
			{
				res[i].add(data[i]);
			}
			else
			{
				res[forest.map[i]].add(data[i]);
			}
		}
		return res;
	}
	
	/*
	 * An edge between two variants indicating that they can be merged
	 * Sorting is non-decreasing order of edge weights (ties broken with variant IDs)
	 */
	class Edge implements Comparable<Edge>
	{
		int from, to;
		double dist;
		Edge(int from, int to, double dist)
		{
			this.from = from;
			this.to = to;
			this.dist = dist;
		}
		@Override
		public int compareTo(Edge o) {
			if(Math.abs(dist - o.dist) > 1e-9)
			{
				return Double.compare(dist, o.dist);
			}
			if(data[from].hash != data[o.from].hash) return data[from].hash - (data[o.from].hash);
			if(data[to].hash != data[o.to].hash) return data[to].hash - (data[o.to].hash);
			if(from != o.from) return data[from].id.compareTo(data[o.from].id);
			return data[to].id.compareTo(data[o.to].id);
		}
	}
}
