/*
 * Main interface for merging variants.  This assumes all variants are on the
 * same chromosome, and have the same type/strand if that separation is desired.
 */
import java.util.ArrayList;
import java.util.PriorityQueue;

public class VariantMerger {
	
	Variant[] data;
	int n;
	double distThreshold;
	Forest forest;
	KDTree knn;

	public VariantMerger(Variant[] data, double distThreshold)
	{
		
		this.distThreshold = distThreshold;
		n = data.length;
		for(int i = 0; i<n; i++) data[i].index = i;
		this.data = data;
		
		forest = new Forest(data);
		knn = new KDTree(data);
	}
		
	void runMerging()
	{
		int[] countEdgesProcessed = new int[n];
		Variant[][] nearestNeighbors = new Variant[n][];
		
		PriorityQueue<Edge> toProcess = new PriorityQueue<Edge>();
		for(int i = 0; i<n; i++)
		{
			nearestNeighbors[i] = knn.kNearestNeighbor(data[i], 4);
			toProcess.add(new Edge(i, nearestNeighbors[i][0].index, data[i].distance(nearestNeighbors[i][0])));
		}
		
		while(!toProcess.isEmpty())
		{
			Edge e = toProcess.poll();
			boolean valid = forest.union(e.from, e.to);
			if(valid)
			{
				System.out.println("Added edge from " + data[e.from].id + " to " + data[e.to].id + " with distance " + e.dist);
			}
			
			while(true)
			{
				countEdgesProcessed[e.from]++;
				if(countEdgesProcessed[e.from]== nearestNeighbors[e.from].length)
				{
					nearestNeighbors[e.from] = knn.kNearestNeighbor(data[e.from], 2 * nearestNeighbors[e.from].length);
				}
				Variant candidateTo = nearestNeighbors[e.from][countEdgesProcessed[e.from]];
				if(data[e.from].distance(candidateTo) > distThreshold + 1e-9)
				{
					// This edge was invalid because of the distance from the query, so stop looking at any edges since they'll only get farther away
					break;
				}
				else if(data[e.from].sample != candidateTo.sample)
				{
					// Possibly valid edge
					toProcess.add(new Edge(e.from, candidateTo.index, data[e.from].distance(candidateTo)));
				}
				else
				{
					// This edge was invalid because of coming from the same sample, so ignore it and try the next one
					continue;
				}
			}
		}
	}
	
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
			if(from != o.from) return data[from].id.compareTo(data[o.from].id);
			return data[to].id.compareTo(data[o.to].id);
		}
	}
}
