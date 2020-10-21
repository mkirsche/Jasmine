import java.io.File;
import java.io.FileInputStream;
import java.util.Scanner;

import java.util.Scanner;

public class test {
public static void main(String[] args) throws Exception
{
	String fn = "/home/mkirsche/Downloads/problematic_line.txt";
	Scanner input = new Scanner(new FileInputStream(new File(fn)));
	String line = input.nextLine();
	System.out.println(line.split("\t").length);
	VcfEntry e = VcfEntry.fromLine(line);
	System.out.println(e.getPos());
}
}
