/*
 * Minimum set of info we need to store for a variant
 */
public class Variant
{
	int sample; // Which sample number the variant came from
	String id; // Variant ID, assumed to be unique for all variants.  It has the "<sampleId>_" added to the beginning to ensure this.
	int start, end; // End should be start+length for insertion
	String graphID; // Store chromosome, and optionally type and strand
	int index; // this is initialized and used internally for bookkeeping and does not come from VCF
	Variant(int sample, String id, int start, int end, String graphID)
	{
		this.sample = sample;
		this.id = id;
		this.start = start;
		this.end = end;
		this.graphID = graphID;
	}
	double distance(Variant v)
	{
		double dStart = start - v.start;
		double dEnd = end - v.end;
		return Math.sqrt(dStart * dStart + dEnd * dEnd); 
	}
	public String toString()
	{
		return "id: " + id + ", sample: " + sample + ", start: " + start + ", end: " + end;
	}
}
