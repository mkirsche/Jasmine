/*
 * Map for normalizing chromosome names
 */
import java.io.File;
import java.io.FileInputStream;
import java.util.HashMap;
import java.util.Scanner;

public class ChrNameNormalization 
{
	HashMap<String, String> normMap = new HashMap<String, String>();
	ChrNameNormalization() throws Exception
	{
		normMap = new HashMap<String, String>();
		if(Settings.DEFAULT_CHR_NORM)
		{
			// Remove "chr" from chromosome names
			for(int i = 1; i<=22; i++)
			{
				normMap.put("chr" + i, i + "");
			}
			normMap.put("chrX", "X");
			normMap.put("chrY", "Y");
			normMap.put("chrM", "MT");
		}
		else if(Settings.CHR_NORM_FILE.length() > 0)
		{
			// Read in chromosome name map
			Scanner input = new Scanner(new FileInputStream(new File(Settings.CHR_NORM_FILE)));
			while(input.hasNext())
			{
				String line = input.nextLine();
				if(line.length() == 0)
				{
					continue;
				}
				String[] tokens = line.split(" ");
				for(int i = 1; i<tokens.length; i++) normMap.put(tokens[i], tokens[0]);
			}
		}
	}
	
	/*
	 * Returns normalized chromosome name
	 */
	String normalize(String chrName)
	{
		if(normMap.containsKey(chrName))
		{
			return normMap.get(chrName);
		}
		return chrName;
	}
}
