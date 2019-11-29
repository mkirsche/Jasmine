/*
 * Class for testing variant merging
*/

import java.util.ArrayList;

public class VariantMergeTest {
	public static void main(String[] args)
	{
		Variant[] data = new Variant[] {
				new Variant(0, "var1", 10, 5, "chr1", null),
				new Variant(0, "var2", 1, 5, "chr1", null),
				new Variant(0, "var3", 18, 5, "chr1", null),
				new Variant(1, "var4", 12, 7, "chr1", null),
				new Variant(1, "var5", 10, 5, "chr1", null),
				new Variant(1, "var6", 30, 30, "chr1", null),
				new Variant(1, "var7", 0, 0, "chr1", null),
				new Variant(2, "var8", 12, 12, "chr1", null),
				new Variant(2, "var9", 15, 15, "chr1", null),
				new Variant(2, "var10", 20, 20, "chr1", null),
				new Variant(2, "var11", 28, 28, "chr1", null),
				new Variant(3, "var12", 25, 25, "chr1", null),
				new Variant(4, "var13", 22, 22, "chr1", null)
		};
		
		Settings.MAX_DIST = 5;
		VariantMerger vm = new VariantMerger(data);
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
