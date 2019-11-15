/*

	Main class for testing variant merging.

	Current output of this (only prints merged variants):
	
	Variant:
	id: var1, sample: 0, start: 10, end: 5
	id: var5, sample: 1, start: 10, end: 5
	
	Variant:
	id: var6, sample: 1, start: 30, end: 30
	id: var11, sample: 2, start: 28, end: 28
	id: var12, sample: 3, start: 25, end: 25
	
	Variant:
	id: var4, sample: 1, start: 12, end: 7
	id: var8, sample: 2, start: 12, end: 12
	
	Variant:
	id: var10, sample: 2, start: 20, end: 20
	id: var13, sample: 4, start: 22, end: 22

*/

import java.util.ArrayList;

public class Thriver {
	public static void main(String[] args)
	{
		Variant[] data = new Variant[] {
				new Variant(0, "var1", 10, 5, "chr1"),
				new Variant(0, "var2", 1, 5, "chr1"),
				new Variant(0, "var3", 18, 5, "chr1"),
				new Variant(1, "var4", 12, 7, "chr1"),
				new Variant(1, "var5", 10, 5, "chr1"),
				new Variant(1, "var6", 30, 30, "chr1"),
				new Variant(1, "var7", 0, 0, "chr1"),
				new Variant(2, "var8", 12, 12, "chr1"),
				new Variant(2, "var9", 15, 15, "chr1"),
				new Variant(2, "var10", 20, 20, "chr1"),
				new Variant(2, "var11", 28, 28, "chr1"),
				new Variant(3, "var12", 25, 25, "chr1"),
				new Variant(4, "var13", 22, 22, "chr1")
		};
		
		VariantMerger vm = new VariantMerger(data, 5.0);
		vm.runMerging();
		ArrayList<Variant>[] res = vm.getGroups();
		
		System.out.println();
		for(ArrayList<Variant> list : res)
		{
			if(list.size() > 1)
			{
				System.out.println("Variant:");
				for(Variant v : list)
				{
					System.out.println(v);
				}
				System.out.println();
			}
		}
	}
}
