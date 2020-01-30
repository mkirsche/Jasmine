/*
 * Multi-threading support for variant merging
 * Since each chromosome (and possibly type and strand) is its own graph,
 * the algorithm can be parallelized pretty naturally.
 * 
 * The graphs that need to be processed are stored in a queue, and each thread
 * processes one graph at a time, querying the queue for the next graph to process
 */

import java.util.ArrayList;
import java.util.Collections;
import java.util.TreeMap;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.atomic.AtomicInteger;

public class ParallelMerger {
	
	// IDs of graphs left to process
	ConcurrentLinkedQueue<String> todo;
	
	// the variant graphs on which merging will be performed
	TreeMap<String, ArrayList<Variant>> allVariants;
	
	// A data structure for holding merged variants to output
	VariantOutput output;
	
	// The number of threads to use
	int numThreads;
	
	// The total number of samples 
	int sampleCount;
	
	AtomicInteger totalMerged = new AtomicInteger(0);
	
	ParallelMerger(TreeMap<String, ArrayList<Variant>> allVariants, VariantOutput output, int sampleCount)
	{
		this.allVariants = allVariants;
		this.output = output;
		this.numThreads = Settings.THREADS;
		System.out.println("Nummber of threads: " + numThreads);
		this.sampleCount = sampleCount;
		todo = new ConcurrentLinkedQueue<String>();
		for(String s : allVariants.keySet())
		{
			todo.add(s);
		}
	}
	
	/*
	 * Start merging in parallel, initializing all threads
	 */
	void run() throws Exception
	{
		// The last thread in the array is the main thread, so it calls
		// run() instead of start() and doesn't get joined below
		MyThread[] threads = new MyThread[numThreads];
		for(int i = 0; i<numThreads; i++)
		{
			threads[i] = new MyThread();
			if(i == numThreads - 1)
			{
				threads[i].run();
			}
			else
			{
				threads[i].start();
			}
		}
		for(int i = 0; i<numThreads-1; i++)
		{
			threads[i].join();
		}
	}

	/*
	 * A single thread performing variant merging
	 */
	public class MyThread extends Thread {
			
		public void run()
		{
			while(!todo.isEmpty())
			{
				String graphID = todo.poll();
				System.out.println("Merging graph ID: " + graphID);
				ArrayList<Variant> variantList = allVariants.get(graphID);
				Collections.sort(variantList);
				VariantMerger vm = new VariantMerger(variantList);
				vm.runMerging();
				ArrayList<Variant>[] res = vm.getGroups();
				output.addGraph(graphID, res, sampleCount);
				int merges = 0;
				for(ArrayList<Variant> list : res)
				{
					if(list.size() > 1)
					{
						merges++;
					}
				}
				totalMerged.addAndGet(merges);
			}
			
		}
		
	}
}
