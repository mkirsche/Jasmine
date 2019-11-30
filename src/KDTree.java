/*
 * Data structure for fast k-nearest neighbor queries in variant sets
 * For a given query, the k closest points to it in the dataset will be reported,
 * breaking ties by variant ID to ensure deterministic behavior.
 * 
 * We assume variants are 2-D points; nearness is based on Euclidean distance or its generalizations.
 * 
 * Uses algorithm described here: 
 * https://courses.cs.washington.edu/courses/cse599c1/13wi/slides/lsh-hashkernels-annotated.pdf
 */

import java.util.Collections;
import java.util.LinkedList;
import java.util.PriorityQueue;

public class KDTree 
{
	Node root;
	Node search;
	PriorityQueue<Candidate> best;
	int cnt;
	int querySize;
	int K;
	
	/*
	 * Initializes a KD-tree from a list of variants
	 */
	public KDTree(Variant[] p) 
	{
		K = 2;
		LinkedList<Node> list = new LinkedList<Node>();
		for (Variant q : p) list.add(new Node(q));
		root = build(list, 0);
	}
	
	/*
	 * Build a subtree of the data structure from a subset of the points
	 * Initially this is called on the whole tree with all the points
	 */
	private Node build(LinkedList<Node> p, int depth) 
	{
		if (p.size() == 0) return null;
		Node pivot = p.remove();
		
		// Sort the points into left and right subtrees based on current split dimension
		LinkedList<Node> left = new LinkedList<Node>();
		LinkedList<Node> right = new LinkedList<Node>();
		while (!p.isEmpty()) 
		{
			if (p.peek().planes[depth % K] < pivot.planes[depth % K])
				left.add(p.remove());
			else
				right.add(p.remove());
		}
		pivot.children[0] = build(left, depth + 1);
		pivot.children[1] = build(right, depth + 1);
		
		return pivot;
	}
	
	/*
	 * Gets the k nearest neighbors for a query variant
	 */
	public Variant[] kNearestNeighbor(Variant p, int k) {
		search = new Node(p);
		best = new PriorityQueue<Candidate>();
		querySize = k;
		search(root, 0);
		Variant[] res = new Variant[best.size()];
		int idx = res.length - 1;
		while(!best.isEmpty())
		{
			res[idx--] = best.poll().v;
		}
		return res;
	}
	
	/*
	 * Search the subtree rooted at cur for candidate points in the set of query's k-nearest neighbors
	 */
	private void search(Node cur, int depth) {
		if (cur == null) return;
		int betterChild = (int) Math.signum(search.planes[depth % K] - cur.planes[depth % K]) < 0 ? 0 : 1;
		search(cur.children[betterChild], depth + 1);
		Candidate toAdd = new Candidate(cur.p, cur.p.distance(search.p));
		if (best == null || best.size() < querySize || toAdd.compareTo(best.peek()) > 0) 
		{
			if(best.size() == querySize)
			{
				best.poll();
			}
			best.add(toAdd);
		}
		if (best.size() < querySize || Math.abs(search.planes[depth % K] - cur.planes[depth % K]) < best.peek().dist)
		{
			search(cur.children[1 - betterChild], depth + 1);
		}
	}
	
	/*
	 * A node of the KD tree
	 * Each node has a variant, storing alongside it its values along the split planes, as well as two (possibly null) children
	 */
	private class Node {
		Node[] children;
		Variant p;
		double[] planes;
		public Node(Variant pp) 
		{
			p = pp;
			planes = new double[K];
			planes[0] = p.start;
			planes[1] = p.end; // add additional dimensions as necessary
			children = new Node[2];
		}
	}
	
	/*
	 * Candidate k-nearest neighbor of the current query point
	 */
	private static class Candidate implements Comparable<Candidate>
	{
		Variant v;
		double dist;
		Candidate(Variant v, double dist)
		{
			this.v = v;
			this.dist = dist;
		}
		public int compareTo(Candidate o)
		{
			if(Math.abs(dist - o.dist) > 1e-9) return Double.compare(o.dist, dist);
			if(v.hash != o.v.hash) return o.v.hash - v.hash;
			return o.v.id.compareTo(v.id);
			
		}
	}
}
