/*
 * Minimum set of info we need to store for a variant
 */

public class Variant implements Comparable<Variant>
{
	// Which sample number the variant came from
	int sample; 
	
	// Variant ID, assumed to be unique for all variants.  It has the "<sampleId>_" added to the beginning to ensure this.
	String id; 
	
	// End should be start+length for insertion
	double start, end;
	
	// Store chromosome, and optionally type and strand
	String graphID;
	
	// This is initialized and used internally for bookkeeping and does not come from VCF
	int index;
	
	// For insertions, the sequence being inserted, or null otherwise
	String seq;
	
	// The maximum distance a variant can be away to merge with this one
	int maxDist;
	
	// The minimum sequence similarity another variant needs to merge with this one if both are insertions
	double minSeqId;
	
	// The stat and end of an interval for checking overlap
	double[] interval;
	
	int hash;
	static int hash(String infoString)
	{
		long res = 0;
		int mod = (int)(1e9+7);
		char[] cs = infoString.toCharArray();
		for(char c : cs)
		{
			res = ((res * 17) + c)%mod;
		}
		return (int)res;
	}
	
	/*
	 * Returns the distance from the variant's (start, end) pair to a given (x, y) point
	 */
	double distFromPoint(double x, double y)
	{
		double dStart = start - x;
		double dEnd = end - y;
		int norm = Settings.KD_TREE_NORM;
		if(norm == 2)
		{
			return Math.sqrt(dStart * dStart + dEnd * dEnd); 
		}
		else
		{
			double powSum = Math.abs(Math.pow(dStart, norm)) + Math.abs(Math.pow(dEnd, norm));
			return Math.pow(powSum, 1.0 / norm);
		}
	}
	
	Variant(int sample, String id, double start, double end, String graphID, String seq, int maxDist, double minSeqId)
	{
		this.sample = sample;
		this.id = id;
		this.start = start;
		this.end = end;
		this.graphID = graphID;
		if(minSeqId > 0) this.seq = seq;
		this.maxDist = maxDist;
		this.minSeqId = minSeqId;
		hash = 0;
		interval = null;
	}
	
	Variant(int sample, String id, double start, double end, String graphID, String seq)
	{
		this.sample = sample;
		this.id = id;
		this.start = start;
		this.end = end;
		this.graphID = graphID;
		if(minSeqId > 0) this.seq = seq;
		this.maxDist = Settings.MAX_DIST;
		this.minSeqId = Settings.MIN_SEQUENCE_SIMILARITY;
		interval = null;
	}
	
	/*
	 * The distance between two variants based on differences in their start/end coordinates
	 * The metric used is a generalization of Euclidean distance
	 */
	double distance(Variant v)
	{
		return distFromPoint(v.start, v.end);
	}
	
	/*
	 * The similarity score of two variants, which is based on sequence similarity for pairs of insertions and 1 otherwise
	 */
	double stringSimilarity(Variant v)
	{
		// If either sequence is null, the variant is either non-insertion, or has no sequence, which we don't want to penalize for
		if(seq == null || v.seq == null)
		{
			return 1.0;
		}
		
		if(Settings.USE_EDIT_DISTANCE)
		{
			return StringUtils.editDistanceSimilarity(seq, v.seq);
		}
		else
		{
			return StringUtils.jaccardSimilarity(seq, v.seq);
		}
	}
	
	/*
	 * Whether or not the sequence similarity of two variants is high enough for them to be merged
	 */
	boolean passesStringSimilarity(Variant v)
	{
		if(seq == null || v.seq == null)
		{
			return true;
		}
		
		String s = seq, t = v.seq;
		
		double similarityNeeded = Math.min(minSeqId, v.minSeqId);
		
		// If there is no sequence identity requirement, don't compute the score and just return true
		if(similarityNeeded <= 0)
		{
			return true;
		}
		
		int minLength = Math.min(s.length(), t.length());
		int maxLength = s.length() + t.length() - minLength;
		
		if(minLength < maxLength * similarityNeeded - 1e-9)
		{
			return false;
		}
		
		return stringSimilarity(v) >= similarityNeeded - 1e-9;
	}
	
	/*
	 * Human-readable format for printing some of the variant information
	 */
	public String toString()
	{
		return "id: " + id + ", sample: " + sample + ", start: " + start + ", end: " + end;
	}

	public int compareTo(Variant o) 
	{
		if(hash != o.hash) return Long.compare(hash, o.hash);
		if(start != o.start) return Double.compare(start, o.start);
		return id.compareTo(o.id);
	}

	public boolean passesOverlap(Variant v) 
	{
		if(interval == null || v.interval == null)
		{
			return true;
		}
		double maxStart = Math.max(interval[0], v.interval[0]);
		double minEnd = Math.min(interval[1], v.interval[1]);
		if(minEnd <= maxStart + 1E-9)
		{
			return false;
		}
		double maxIntervalSize = Math.max(interval[1] - interval[0], v.interval[1] - v.interval[0]);
		return minEnd - maxStart + 1e-9 >= maxIntervalSize * Settings.OVERLAP_REQUIRED;
	}
}
