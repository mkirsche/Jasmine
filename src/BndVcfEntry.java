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
		return Long.parseLong(altTokens[1].split(":")[1]);
	}
	
	/*
	 * Always call the type of a BND-style VCF translocation to ensure consistently with different formats
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
		return getTranslocationGraphID();
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
		return altTokens[1].split(":")[0];
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
		return "";
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
