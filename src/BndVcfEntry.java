/*
 * A VCF entry where the ALT field contains most of the SV information in BND format
 * It is assumed that all such SVs are translocations and that each entry encodes its own variant
 */
public class BndVcfEntry extends VcfEntry {
	
	String[] altTokens;
	public BndVcfEntry(String line) throws Exception
	{
		super(line);
		altTokens =  getAlt().split("[\\[\\]]");
	}
	
	/*
	 * Fall back on the length being zero if the SVLEN field is not there
	 */
	public int getLength() throws Exception
	{
		if(hasInfoField("SVLEN"))
		{
			return Integer.parseInt(getInfo("SVLEN"));
		}
		return 0;
	}
	
	/*
	 * Getting the end coordinate may require parsing the ALT field
	 */
	public long getEnd() throws Exception
	{
		if(hasInfoField("END"))
		{
			return Long.parseLong(getInfo("END"));
		}
		String chrPosToken = altTokens[1];
		return Long.parseLong(chrPosToken.substring(1 + chrPosToken.lastIndexOf(':')));
	}
	
	/*
	 * Always call the type of a BND-style VCF translocation to ensure consistency with different formats
	 */
	public String getType()
	{
		return "TRA";
	}
	
	/*
	 * The graph ID here has to consider both chromosomes
	 */
	public String getGraphId() throws Exception
	{
		String first = getChromosome(), second = getChr2();
		if(first.compareTo(second) > 0)
		{
			String tmp = first;
			first = second;
			second = tmp;
		}
		String id = first + "_" + second;
		if(Settings.USE_TYPE)
		{
			id += "_" + getType();
		}
		if(Settings.USE_STRAND)
		{
			id += "_" + getStrand();
		}
		return id;
	}
	
	/*
	 * The second chromosome can be found in either the CHR2 INFO field or the ALT field
	 */
	public String getChr2() throws Exception
	{
		if(hasInfoField("CHR2"))
		{
			return getInfo("CHR2");
		}
		String chrPosToken = altTokens[1];
		return chrPosToken.substring(0, chrPosToken.lastIndexOf(':'));
	}
	
	/*
	 * The strands may need to be inferred from the ALT square bracket format
	 */
	public String getStrand() throws Exception
	{
		String res = getInfo("STRANDS");
		if(res.length() == 0)
		{
			return strandsFromAltFormat();
		}
		return res;
	}
	
	/*
	 * Determine the strands from the ALT square bracket format
	 */
	public String strandsFromAltFormat()
	{
		String alt = getAlt();
		if(alt.startsWith("["))
		{
			return "+-";
		}
		else if(alt.startsWith("]"))
		{
			return "--";
		}
		else if(alt.contains("["))
		{
			return "++";
		}
		else if(alt.contains("]"))
		{
			return "-+";
		}
		return "";
	}
	
	/*
	 * Gets the first coordinate of the variant
	 */
	public double getFirstCoord() throws Exception
	{
		if(getChromosome().compareTo(getChr2()) > 0)
		{
			if(hasInfoField("AVG_END"))
			{
				return Double.parseDouble(getInfo("AVG_END"));
			}
			return getEnd();
		}
		if(hasInfoField("AVG_START"))
		{
			return Double.parseDouble(getInfo("AVG_START"));
		}
		return getPos();
	}
	
	/*
	 * Since length is undefined, get the second coord instead
	 */
	public double getSecondCoord() throws Exception
	{
		if(getChromosome().compareTo(getChr2()) > 0)
		{
			if(hasInfoField("AVG_START"))
			{
				return Double.parseDouble(getInfo("AVG_START"));
			}
			return getPos();
		}
		if(hasInfoField("AVG_END"))
		{
			return Double.parseDouble(getInfo("AVG_END"));
		}
		return getEnd();
	}
}
