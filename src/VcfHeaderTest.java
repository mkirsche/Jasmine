/*
 * Test to make sure VCF header is properly handling adding INFO fields and adding/resetting FORMAT fields
 * Output should be (in order): Test, SVTYPE, INFO1-3, FORMAT1-3, and #CHR
 */
import java.io.PrintWriter;

public class VcfHeaderTest {
public static void main(String[] args)
{
	VcfHeader header = new VcfHeader();
	header.addLine("#Test");
	header.addLine("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of the SV.\">");
	header.addLine("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
	header.addLine("##CHR etc.");
	header.addInfoField("INFO1", "1", "String", "desc1");
	header.addInfoField("INFO2", "1", "String", "desc2");
	header.addFormatField("FORMAT3", "1", "String", "descf3");
	header.resetFormatFields();
	header.addFormatField("FORMAT1", "1", "String", "descf1");
	header.addInfoField("INFO3", "1", "String", "desc3");
	header.addFormatField("FORMAT2", "1", "String", "descf2");
	header.addFormatField("FORMAT3", "1", "String", "descf3");
	header.addFormatField("FORMAT1", "1", "String", "descf1");
	PrintWriter out = new PrintWriter(System.out);
	header.print(out);
	out.close();
}
}
