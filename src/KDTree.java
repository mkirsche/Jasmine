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

import java.util.ArrayDeque;
import java.util.ArrayList;
import java.util.Arrays;
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
	
	int n;
	
	/*
	 * Initializes a KD-tree from a list of variants
	 */
	public KDTree(Variant[] p) 
	{
		n = p.length;
		K = 2;
		LinkedList<Node> list = new LinkedList<Node>();
		for (Variant q : p) list.add(new Node(q));
		root = build(list, 0);//buildNonrecursive(list, 0).get(0);
	}
	
	public KDTree(Variant[] p, boolean recursive) 
	{
		n = p.length;
		K = 2;
		LinkedList<Node> list = new LinkedList<Node>();
		for (Variant q : p) list.add(new Node(q));
		root = recursive ? build(list, 0) : buildNonrecursive(list).get(0);
	}
	
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
	 * Build the data structure from a list of points without recursion
	 * This avoids stack overflow issues caused by larger datasets
	 */
	private ArrayList<Node> buildNonrecursive(LinkedList<Node> p) 
	{
		ArrayList<Node> nodeList = new ArrayList<Node>();
		if(p.size() == 0)
		{
			return null;
		}
		
		// The stack of node lists to process (in place of recursive calls)
		ArrayDeque<LinkedList<Node>> toProcess = new ArrayDeque<LinkedList<Node>>();
		ArrayDeque<Integer> parents = new ArrayDeque<Integer>();
		ArrayDeque<Integer> depths = new ArrayDeque<Integer>();
		ArrayDeque<Integer> parentsides = new ArrayDeque<Integer>();
		
		// Initialize root to null to be filled 
		//nodeList.add(res);
		toProcess.addFirst(p);
		parents.addFirst(-1); // This is not actually the parent of the root, but it will be ignored anyways
		depths.addFirst(0);
		parentsides.addFirst(-1);
		
		while(!toProcess.isEmpty())
		{
			// Get information for processing this node from stacks
			LinkedList<Node> pcur = toProcess.pollFirst();
			int parentcur = parents.pollFirst();
			int depthcur = depths.pollFirst();
			int parentsidecur = parentsides.pollFirst();
			
			// Get pivot as the first point in the list
			Node pivot = pcur.remove();
			
			// Separate this point into points left of the pivot vs. right of the pivot
			LinkedList<Node> left = new LinkedList<Node>();
			LinkedList<Node> right = new LinkedList<Node>();
			while (!pcur.isEmpty()) 
			{
				Node check = pcur.pollFirst();
				if (check.planes[depthcur % K] < pivot.planes[depthcur % K])
					left.add(check);
				else
					right.add(check);
			}
			
			//pcur.clear();
			
			// Update this node's parent's child-pointer to this node.
			if(parentsidecur != -1)
			{
				nodeList.get(parentcur).children[parentsidecur] = pivot;
			}
			
			pivot.children[0] = null;
			pivot.children[1] = null;
			nodeList.add(pivot);
			
			// Add right child to processing stack
			if(right.size() > 0)
			{
				toProcess.addFirst(right);
				parents.addFirst(nodeList.size() - 1);
				parentsides.addFirst(1);
				depths.addFirst(depthcur + 1);
			}
			
			// Add left child to processing stack
			if(left.size() > 0)
			{
				toProcess.addFirst(left);
				parents.addFirst(nodeList.size() - 1);
				parentsides.addFirst(0);
				depths.addFirst(depthcur + 1);
			}
		}
		
		return nodeList;
	}
	
	/*
	 * Used to make sure two KD-trees are the same
	 */
	static boolean compare(String pref, Node a, Node b)
	{
		if(a == null && b != null)
		{
			System.out.println(pref + " only a is null");
			return true;
		}
		if(b == null && a != null)
		{
			System.out.println(pref + " only b is null");
			return true;
		}
		if(a == null && b == null)
		{
			return false;
		}
		if(a.planes[0] != b.planes[0] || a.planes[1] != b.planes[1])
		{
			System.out.println(pref + " diff value: " + a.planes[0] + " " + a.planes[1] + " " + b.planes[0] + " " + b.planes[1]);
			return true;
		}
		
		boolean leftDiff = compare(pref + "L", a.children[0], b.children[0]);
		if(leftDiff)
		{
			return true;
		}
		else
		{
			return compare(pref + "R", a.children[1], b.children[1]);
		}
		
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
			children[0] = null;
			children[1] = null;
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
