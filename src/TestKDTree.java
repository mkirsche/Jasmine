/*
 * Basic test for making sure my KD tree is behaving reasonably
 */

public class TestKDTree {
public static void main(String[] args)
{
	Variant[] data = new Variant[] {
			new Variant(0, "var1", 10, 5, "chr1", null),
			new Variant(0, "var2", 1, 5, "chr1", null),
			new Variant(0, "var3", 18, 5, "chr1", null),
			new Variant(0, "var4", 12, 7, "chr1", null),
			new Variant(0, "var5", 10, 5, "chr1", null),
			new Variant(0, "var6", 30, 30, "chr1", null),
			new Variant(0, "var7", 0, 0, "chr1", null)
	};
	
	KDTree kdt = new KDTree(data);
	for(int i = 0; i<data.length;i++)
	{
		Variant[] cur = kdt.kNearestNeighbor(data[i], 5);
		for(Variant v : cur)
		{
			System.out.println(v.id+" "+v.start+" "+v.end+" "+v.distance(data[i]));
		}
		System.out.println();
	}
}
}
