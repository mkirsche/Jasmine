/*
 * Methods for storing VCF v4.2 entries for structural variants
 * It allows easy access to main variant fields, plus some methods
 * to do more parsing and error-checking for things like type and length.
 */

import java.util.Arrays;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class VcfEntry {

	String originalLine;
	String[] tabTokens;
	String oldId;
	
	public static VcfEntry fromLine(String line) throws Exception
	{
		VcfEntry res = new VcfEntry(line);
		if(res.getType().equals("BND"))
		{
			return new BndVcfEntry(line);
		}
		
		// Adding this for reverse compatibility
		if(res.hasInfoField("IN_SPECIFIC"))
		{
			res.setInfo("IS_SPECIFIC", res.getInfo("IN_SPECIFIC"));
		}
		return res;
	}
	
	public VcfEntry(String line) throws Exception
	{
		originalLine = line;
		tabTokens = line.split("\t");
		if(tabTokens.length < 8)
		{
			throw new Exception("VCF line had too few entries: "
					+ Arrays.toString(tabTokens));
		}
		oldId = getId();
		if(Settings.NORMALIZE_TYPE && !getType().equals("BND"))
		{
			setType(getNormalizedType());
		}
		setChromosome(Settings.CHR_NAME_MAP.normalize(getChromosome()));
	}
	
	/*
	 * Update the breakpoint and length fields to be consistent
	 */
	void fixImprecision() throws Exception
	{
		String normalizedType = getNormalizedType();
		
		if(!normalizedType.equals("TRA") && !getAlt().contains("<"))
		{
			setInfo("SVLEN", getAlt().length() - getRef().length() + "");
		}
		
		if(normalizedType.equals("INS"))
		{
			setInfo("END", getPos()+"");
		}
		else if(normalizedType.equals("DEL") || normalizedType.equals("INV") || normalizedType.equals("DUP"))
		{
			// If no indicator of length, try using the end to fill it
			if(getAlt().contains("<") && !hasInfoField("SVLEN") && !hasInfoField("SEQ"))
			{
				if(hasInfoField("END"))
				{
					if(normalizedType.equals("DEL"))
					{
						setInfo("SVLEN", getPos() - Long.parseLong(getInfo("END")) + "");
					}
					else
					{
						setInfo("SVLEN", Long.parseLong(getInfo("END")) - getPos() + "");
					}
				}
			}
			
			// If there is a length indicator, use it to set/fix end if necessary
			else
			{
				setInfo("END", getPos() + Math.abs(getLength())+"");
			}
		}
	}
	
	/*
	 * Reconstruct the VCF line by concatenating and tab-separating the fields
	 */
	public String toString()
	{
		StringBuilder sb = new StringBuilder("");
		for(int i = 0; i<tabTokens.length; i++)
		{
			sb.append(tabTokens[i]);
			if(i < tabTokens.length - 1)
			{
				sb.append("\t");
			}
		}
		return sb.toString();
	}
	
	/*
	 * Get the chromosome field
	 */
	public String getChromosome()
	{
		return tabTokens[0];
	}
	
	/*
	 * Set the chromosome field
	 */
	public void setChromosome(String s)
	{
		tabTokens[0] = s;
	}
	
	/*
	 * Get the POS field
	 */
	public long getAvgPos() throws Exception
	{
		if(hasInfoField("AVG_START"))
		{
			return (long)(.5 + Double.parseDouble(getInfo("AVG_START")));
		}
		return getPos();
	}
	
	/*
	 * Get the POS field
	 */
	public long getPos() throws Exception
	{
		return Long.parseLong(tabTokens[1]);
	}
	
	/*
	 * Set the POS field
	 */
	public void setPos(long val)
	{
		tabTokens[1] = val+"";
	}
	
	/*
	 * Get the variant ID field
	 */
	public String getId()
	{
		return tabTokens[2];
	}
	
	/*
	 * Set the variant ID field
	 */
	public void setId(String s)
	{
		tabTokens[2] = s;
	}
	
	/*
	 * Get the REF sequence field
	 */
	public String getRef()
	{
		return tabTokens[3];
	}
	
	/*
	 * Set the REF sequence field
	 */
	public void setRef(String s)
	{
		tabTokens[3] = s;
	}
	
	/*
	 * Get the ALT sequence field
	 */
	public String getAlt()
	{
		return tabTokens[4];
	}
	
	/*
	 * Set the ALT sequence field
	 */
	public void setAlt(String s)
	{
		tabTokens[4] = s;
	}
	
	/*
	 * The length of a variant
	 */
	public int getLength() throws Exception
	{
		try {
			String s = getInfo("SVLEN");
			double val = Double.parseDouble(s);
			if(val < 0)
			{
				return (int)(val - .5);
			}
			else
			{
				return (int)(val + .5);
			}
		} catch(Exception e) {
            String seq = getSeq();
            String type = getType();
            if(type.equals("INS")) return seq.length();
            else return -seq.length();
		}
	}
	
	/*
	 * The end position of a variant
	 */
	public long getEnd() throws Exception
	{
		if(hasInfoField("END")) return Long.parseLong(getInfo("END"));
		String type = getType();
		if(type.equals("INS"))
		{
			return getPos();
		}
		else
		{
			return getPos() + Math.abs(getLength());
		}
	}
	
	/*
	 * The end position of a variant
	 */
	public long getAvgEnd() throws Exception
	{
		if(hasInfoField("AVG_END"))
		{
			return (long)(.5 + Double.parseDouble(getInfo("AVG_END")));
		}
		return getEnd();
	}
	
	/*
	 * Get the variant type
	 */
	public String getType() throws Exception
	{
		String res = getInfo("SVTYPE");
		if(res.length() == 0)
		{
			int refLength = getRef().length();
			int altLength = getAlt().length();
			if(getAlt().startsWith("<") && getAlt().endsWith(">"))
			{
				return getAlt().substring(1, getAlt().length() - 1);
			}
			if(refLength > altLength)
			{
				return "INS";
			}
			else if(refLength < altLength)
			{
				return "DEL";
			}
			else
			{
				return "";
			}
		}
		else return res;
	}
	
	/*
	 * Set the SVTYPE INFO field to a certain value
	 */
	public void setType(String s) throws Exception
	{
		setInfo("SVTYPE", s);
	}
	
	/*
	 * The strands on which the variant occurs
	 */
	public String getStrand() throws Exception
	{
		return getInfo("STRANDS");
	}
	
	/*
	 * Get which graph the variant belongs in based on how the user wants them separated
	 */
	public String getGraphID() throws Exception
	{
		if(getType().equals("TRA"))
		{
			return getTranslocationGraphID();
		}
		String id = getChromosome();
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
	 * Gets the graph ID for a translocation, which relies on the chr2 field as well
	 */
	public String getTranslocationGraphID() throws Exception
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
	 * Gets the second chromosome of a translocation from the CHR2 INFO field
	 */
	public String getChr2() throws Exception
	{
		if(hasInfoField("CHR2"))
		{
			return Settings.CHR_NAME_MAP.normalize(getInfo("CHR2"));
		}
		String alt = getAlt();
		if(alt.contains("[") || alt.contains("]"))
		{
			return Settings.CHR_NAME_MAP.normalize(alt.split("[\\[\\]]")[1].split(":")[0]);
		}
		return "";
	}
	
	/*
	 * The genomic sequence being affected/added/deleted by the variant
	 */
	public String getSeq() throws Exception
	{
		if(hasInfoField("SEQ"))
		{
			return getInfo("SEQ");
		}
		String ref = getRef(), alt = getAlt();
		if(alt.startsWith("<"))
		{
			return "";
		}
		String type = getType();
		
		if(type.equals("BND") || type.equals("TRA"))
		{
			return "";
		}
		
		// If the SV is a deletion, swap REF and ALT and treat as an insertion
		if(type.equals("DEL"))
		{
			String tmp = ref;
			ref = alt;
			alt = tmp;
			type = "INS";
		}
		if(ref.equals("X") || ref.equals("N"))
		{
			return alt;
		}
		else if(alt.equals("X") || ref.equals("N"))
		{
			return ref;
		}
		else if(type.equals("INS"))
		{
			int startPad = 0, endPad = 0;
			int totalPad = ref.length();
			while(startPad + endPad < totalPad && startPad < alt.length())
			{
				if(ref.charAt(startPad) == alt.charAt(startPad))
				{
					startPad++;
				}
				else if(ref.charAt(ref.length() - 1 - endPad) == alt.charAt(alt.length() - 1 - endPad))
				{
					endPad++;
				}
				else
				{
					break;
				}
			}
			if(startPad + endPad == totalPad)
			{
				return alt.substring(startPad, alt.length() - endPad);
			}
			else
			{
				return alt;
			}
		}
		else
		{
			return alt;
		}
	}
	
	/*
	 * Set a particular VCF INFO field, adding the field if it doesn't already exist
	 */
	public void setInfo(String field, String val) throws Exception
	{
		String[] infoFields = tabTokens[7].split(";");
		for(String semitoken : infoFields)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(field))
			{
				String updatedToken = key + "=" + val;
				
				// Special case if this is the first INFO field
				if(tabTokens[7].startsWith(semitoken))
				{
					tabTokens[7] = tabTokens[7].replaceFirst(Pattern.quote(semitoken), Matcher.quoteReplacement(updatedToken));
				}
				else
				{
					tabTokens[7] = tabTokens[7].replace(";" + semitoken, ";" + updatedToken);
				}
				return;
			}
		}
		
		// Field not found, so add it!
		tabTokens[7] += ";" + field + "=" + val;
	}
	
	/*
	 * Get the number of supporting reads
	 */
	public int getReadSupport() throws Exception
	{
		if(hasInfoField("RE"))
		{
			return Integer.parseInt(getInfo("RE"));
		}
		String[] reads = getRnames();
		if(reads.length > 0)
		{
			return reads.length;
		}
		return 0;
	}
	
	/*
	 * Get list of supporting read names - first check for RNAMES field, and then anything containing RNAMES
	 */
	public String[] getRnames() throws Exception
	{
		if(hasInfoField("RNAMES"))
		{
			return getInfo("RNAMES").split(",");
		}
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.toUpperCase().contains("RNAMES"))
			{
				return semitoken.substring(1 + equalIndex).split(",");
			}
		}
		return new String[] {};
	}
	
	/*
	 * Get the value of a particular INFO field
	 */
	public String getInfo(String field) throws Exception
	{
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(field))
			{
				return semitoken.substring(1 + equalIndex);
			}
		}
		return "";
	}
	
	/*
	 * Whether or not the line has a particular INFO field
	 */
	public boolean hasInfoField(String fieldName)
	{
		String infoToken = tabTokens[7];
		String[] semicolonSplit = infoToken.split(";");
		for(String semitoken : semicolonSplit)
		{
			int equalIndex = semitoken.indexOf('=');
			if(equalIndex == -1)
			{
				continue;
			}
			String key = semitoken.substring(0, equalIndex);
			if(key.equals(fieldName))
			{
				return true;
			}
		}
		return false;
	}
	
	/*
	 * Get the first value to use for merging, which is the start in all cases but translocations with chr1 > chr2.
	 */
	public double getFirstCoord() throws Exception
	{
		if(getType().equals("TRA") && getChromosome().compareTo(getChr2()) > 0)
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
		return getAvgPos();
	}
	
	/*
	 * Get the second value to use for merging, which depends on the settings
	 */
	public double getSecondCoord() throws Exception
	{
		if(Settings.USE_END || getType().equals("TRA"))
		{
			if(getType().equals("TRA"))
			{
				
			}
			if(hasInfoField("AVG_END"))
			{
				return Double.parseDouble(getInfo("AVG_END"));
			}
			return (int)getEnd();
		}
		
		else
		{
			if(hasInfoField("AVG_LEN"))
			{
				return Double.parseDouble(getInfo("AVG_LEN"));
			}
			return Math.abs(getLength());
		}
	}
	
	/*
	 * Get one of five types: INS, DEL, DUP, INV, or TRA (or "" if none of them fit)
	 */
	public String getNormalizedType() throws Exception
	{
		String type = getType();
		if(type.equals("TRA") || type.equals("BND"))
		{
			return "TRA";
		}
		if(type.equals("INS") || type.equals("DEL") || type.equals("DUP") || type.equals("INV"))
		{
			return type;
		}
		if(hasInfoField("CHR2") && !getInfo("CHR2").equals(getChromosome()))
		{
			return "TRA";
		}
		if(hasInfoField("STRANDS"))
		{
			String strand = getInfo("STRANDS");
			if(strand.equals("++") || strand.equals("--"))
			{
				return "INV";
			}
			else if(strand.equals("-+"))
			{
				return "DUP";
			}
			else if(strand.equals("+-"))
			{
				int length = getLength();
				if(length > 0) return "INS";
				else if(length < 0) return "DEL";
				else return "";
			}
			else return "";
		}
		return "";
	}
	
}
