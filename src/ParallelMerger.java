import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Queue;
import java.util.TreeMap;
import java.util.concurrent.atomic.AtomicInteger;

public class ParallelMerger {
	
	Queue<String> todo;
	
	TreeMap<String, ArrayList<Variant>> allVariants;
	VariantOutput output;
	
	int numThreads;
	int sampleCount;
	
	AtomicInteger totalMerged = new AtomicInteger(0);
	
	ParallelMerger(TreeMap<String, ArrayList<Variant>> allVariants, VariantOutput output, int sampleCount)
	{
		this.allVariants = allVariants;
		this.output = output;
		this.numThreads = Settings.THREADS;
		this.sampleCount = sampleCount;
		todo = new LinkedList<String>();
		for(String s : allVariants.keySet())
		{
			todo.add(s);
		}
	}
	
	void run() throws Exception
	{
		// Here the last thread in the array is the main thread, so it calls
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

	public class MyThread extends Thread {
			
		@Override
		public void run() {
			while(!todo.isEmpty())
			{
				String graphID = todo.poll();
				ArrayList<Variant> variantList = allVariants.get(graphID);
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
