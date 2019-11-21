import java.util.Arrays;

/*
 * Methods for handling VCF v4.2 entries for structural variants
 */
public class VcfEntry implements Comparable<VcfEntry> {

	String originalLine;
	String[] tabTokens;
	
	public VcfEntry(String line) throws Exception
	{
		originalLine = line;
		tabTokens = line.split("\t");
		if(tabTokens.length < 8)
		{
			throw new Exception("VCF line had too few entries: "
					+ Arrays.toString(tabTokens));
		}
	}
	
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
	
	public String getGenotype()
	{
		if(tabTokens.length < 10)
		{
			return "";
		}
		String res = tabTokens[9];
		if(res.indexOf(':') != -1)
		{
			res = res.substring(0, res.indexOf(':'));
		}
		return res;
	}
	
	public String getChromosome()
	{
		return tabTokens[0];
	}
	
	public void setChromosome(String s)
	{
		tabTokens[0] = s;
	}
	
	public long getPos() throws Exception
	{
		try {
			return Long.parseLong(tabTokens[1]);
		} catch(Exception e) {
			throw new Exception("Tried to access invalid VCF position: " + tabTokens[1]);
		}
	}
	
	public void setPos(long val)
	{
		tabTokens[1] = val+"";
	}
	
	public String getId()
	{
		return tabTokens[2];
	}
	
	public String getRef()
	{
		return tabTokens[3];
	}
	
	public void setRef(String s)
	{
		tabTokens[3] = s;
	}
	
	public String getAlt()
	{
		return tabTokens[4];
	}
	
	public void setAlt(String s)
	{
		tabTokens[4] = s;
	}
	
	public int getLength() throws Exception
	{
		String s = getInfo("SVLEN");
		try {
			return (int)(.5 + Double.parseDouble(s));
		} catch(Exception e) {
            try {
                String seq = getSeq();
                String type = getType();
                if(type.equals("INS")) return seq.length();
                else return -seq.length();
            } catch(Exception f) {
			    throw new Exception("SVLEN field is not an integer: " + s);
            }
		}
	}
	
	public void setLength(int len) throws Exception
	{
		setInfo("SVLEN", len+"");
	}
	
	public String getType() throws Exception
	{
		return getInfo("SVTYPE");
	}
	
	public void setType(String s) throws Exception
	{
		setInfo("SVTYPE", s);
	}
	
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
			while(startPad + endPad < totalPad)
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
		throw new Exception("Trying to access VCF INFO field which is not found: " + field);
	}
	
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
	
	public void setInfo(String field, String val) throws Exception
	{
		try {
			String oldVal = getInfo(field);
			String toReplace = field + '=' + oldVal;
			String replacement = field + '=' + val;
			tabTokens[7] = tabTokens[7].replace(toReplace, replacement);
		} catch(Exception e) {
			tabTokens[7] += ";" + field + "=" + val;
		}
	}
		
	public String getKey() throws Exception
	{
		return getChromosome() + ":" + getPos() + ":" + getType() + ":" + getId();
	}
	
	static String getChrFromKey(String key)
	{
		String[] tokens = key.split(":");
		return tokens[0];
	}
	
	static long getPosFromKey(String key)
	{
		String[] tokens = key.split(":");
		return Long.parseLong(tokens[1]);
	}
	
	static String getTypeFromKey(String key)
	{
		String[] tokens = key.split(":");
		return tokens[2];
	}
	
	@Override
	public int compareTo(VcfEntry o) {
		for(int i = 0; i<tabTokens.length && i < o.tabTokens.length; i++)
		{
			if(!tabTokens[i].equals(o.tabTokens[i]))
			{
				return tabTokens[i].compareTo(o.tabTokens[i]);
			}
		}
		return tabTokens.length - o.tabTokens.length;
	}

}
