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
	String seq; // For insertions, the sequence being inserted, or null otherwise
	int maxDist; // The maximum distance a variant can be to merge with this one
	double minSeqId; // The minimum sequence similarity another variant needs to merge with this one if both are insertions
	Variant(int sample, String id, int start, int end, String graphID, String seq, int maxDist, double minSeqId)
	{
		this.sample = sample;
		this.id = id;
		this.start = start;
		this.end = end;
		this.graphID = graphID;
		this.seq = seq;
		this.maxDist = maxDist;
		this.minSeqId = minSeqId;
	}
	Variant(int sample, String id, int start, int end, String graphID, String seq)
	{
		this.sample = sample;
		this.id = id;
		this.start = start;
		this.end = end;
		this.graphID = graphID;
		this.seq = seq;
		this.maxDist = Settings.MAX_DIST;
		this.minSeqId = Settings.MIN_SEQUENCE_SIMILARITY;
	}
	double distance(Variant v)
	{
		double dStart = start - v.start;
		double dEnd = end - v.end;
		return Math.sqrt(dStart * dStart + dEnd * dEnd); 
	}
	double stringSimilarity(Variant v)
	{
		// If either sequence is null, the variant is either non-insertion, or has no sequence, which we don't want to penalize for
		if(seq == null || v.seq == null)
		{
			return 1.0;
		}
		return StringUtils.jaccardSimilarity(seq, v.seq);
		//return StringUtils.editDistanceSimilarity(seq, v.seq);
	}
	boolean passesStringSimilarity(Variant v)
	{
		if(seq == null || v.seq == null)
		{
			return true;
		}
		String s = seq, t = v.seq;
		
		double similarityNeeded = Math.min(minSeqId, v.minSeqId); 
		
		int minLength = Math.min(s.length(), t.length());
		int maxLength = s.length() + t.length() - minLength;
		if(minLength < maxLength * similarityNeeded - 1e-9)
		{
			return false;
		}
		return stringSimilarity(v) >= similarityNeeded - 1e-9;
	}
	public String toString()
	{
		return "id: " + id + ", sample: " + sample + ", start: " + start + ", end: " + end;
	}
}
