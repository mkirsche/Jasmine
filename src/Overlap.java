/*
 * A program for filtering variants based on their overlap with a list of regions.
 */
import java.io.File;
import java.io.FileInputStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.TreeSet;

public class Overlap
{
	static String vcfFn = "";
	static String bedFn = "";
	static String ofn = "";
	static String FILTER_MODE = "CONTAINED_IN_REGION";
	static String REPORT_MODE = "REMOVE";
	static String reportInfo = "";
	
	static ChrNameNormalization chrNorm;
	static void parseArgs(String[] args)
	{
		for(String arg : args)
		{
			int equalsIdx = arg.indexOf('=');
			if(equalsIdx == -1)
			{
				
			}
			else
			{
				String key = arg.substring(0, equalsIdx);
				String val = arg.substring(1 + equalsIdx);
				if(key.equalsIgnoreCase("vcf_file"))
				{
					vcfFn = val;
				}
				else if(key.equalsIgnoreCase("bed_file"))
				{
					bedFn = val;
				}
				else if(key.equalsIgnoreCase("out_file"))
				{
					ofn = val;
				}
				else if(key.equalsIgnoreCase("info_report"))
				{
					reportInfo = val;
				}
			}
		}
		
		if(reportInfo.length() > 0)
		{
			REPORT_MODE = "INFO";
		}
		
		if(vcfFn.length() == 0 || bedFn.length() == 0 || ofn.length() == 0)
		{
			usage();
			System.exit(0);
		}
	}
	static void usage()
	{
		System.out.println();
		System.out.println("Jasmine IGV Screenshot Maker");
		System.out.println("Usage: overlap_jasmine [args]");
		System.out.println("  Example: overlap_jasmine vcf_file=merged.vcf bed_file=regions.bed out_fie=filtered.vcf");
		System.out.println();
		System.out.println("Required args:");
		System.out.println("  vcf_file   (String) - the VCF file with merged SVs");
		System.out.println("  bed_file   (String) - a BED file with regions of interest");
		System.out.println("  out_file   (String) - the name of the output VCF filtered by regions of interest");
		System.out.println();
		System.out.println("Optional args:");
		System.out.println("  info_report (String) [] - the INFO field to indicate presence in regions instead of removing non-overlapping variants");
		System.out.println();
	}
	
	public static void main(String[] args) throws Exception
	{
		parseArgs(args);
		
		Settings.DEFAULT_CHR_NORM = true;
		chrNorm = new ChrNameNormalization();
		
		filterVcf();
	}
	
	static ArrayList<Event> getBedEvents() throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(bedFn)));
		
		ArrayList<Event> events = new ArrayList<Event>();
		int idNum = 0;
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith("#"))
			{
				continue;
			}
			String[] tokens = line.split("\t");
			String chr = tokens[0];
			chr = chrNorm.normalize(chr);
			int start = Integer.parseInt(tokens[1]);
			int end = Integer.parseInt(tokens[2]);
			idNum++;
			events.add(new Event(chr, start, 1, idNum + ""));
			events.add(new Event(chr, end, -1, idNum + ""));
		}
		input.close();
		
		Collections.sort(events);
		
		return events;
	}
	
	/*
	 * Gets the start and end events for variants
	 * Translocations are a special case and are broken up into two event pairs, one for each breakpoint
	 */
	static ArrayList<Event> getVcfEvents() throws Exception
	{
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		
		ArrayList<Event> events = new ArrayList<Event>();
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.startsWith("#"))
			{
				continue;
			}
			VcfEntry entry = VcfEntry.fromLine(line);
			String chr = entry.getChromosome();
			chr = chrNorm.normalize(chr);
			int start = (int)(entry.getPos());
			int end = (int)(entry.getEnd());
			String id = entry.getId();
			
			if(entry.getNormalizedType().equals("TRA"))
			{
				events.add(new Event(chr, start, 1, id + "_breakpoint1"));
				events.add(new Event(chr, start+1, -1, id + "_breakpoint1"));
				String chr2 = entry.getChr2();
				if(chr2.length() != 0)
				{
					events.add(new Event(chr2, end, 1, id + "_breakpoint2"));
					events.add(new Event(chr2, end + 1, -1, id + "_breakpoint2"));
				}
			}
			else
			{
				events.add(new Event(chr, start, 1, id));
				events.add(new Event(chr, end + 1, -1, id));
			}
		}
		input.close();
		
		Collections.sort(events);
		
		return events;
	}
	
	/*
	 * Gets a list of overlaps based on variant and region start and end events
	 * For each variant, it outputs a list of the IDs of regions with which it overlaps
	 */
	static HashMap<String, HashSet<String>> getOverlaps(ArrayList<Event> regions, ArrayList<Event> variants)
	{
		// The next region and variant events to consider
		int regionIdx = 0, variantIdx = 0;
		
		// As we do the plane sweep, the list of regions we are currently inside, if any
		TreeSet<String> openRegions = new TreeSet<String>();
		// The list of variant intervals we are currently inside
		TreeSet<String> openVariants = new TreeSet<String>();

		HashMap<String, HashSet<String>> overlaps = new HashMap<String, HashSet<String>>();
		
		// Plane sweep time!
		while(true)
		{
			// Stop when there are no more events
			if(regionIdx == regions.size() && variantIdx == variants.size())
			{
				break;
			}
			
			// Whether or not the next event to process is a region event (as opposed to a variant)
			boolean nextEventRegion = false;
			
			// If we are out of regions, take a variant event
			if(regionIdx == regions.size())
			{
				nextEventRegion = false;
			}
			
			// If we are out of variants, take a region event
			else if(variantIdx == variants.size())
			{
				nextEventRegion = true;
			}
			
			// If we have both left, compare positions and choose according to overlap scheme
			else
			{
				Event nextRegion = regions.get(regionIdx);
				Event nextVariant = variants.get(variantIdx);
				
				// If region is on earlier chromosome, take that
				if(nextRegion.chr.compareTo(nextVariant.chr) < 0)
				{
					nextEventRegion = true;
				}
				
				// If region is on later chromosome, take variant
				else if(nextRegion.chr.compareTo(nextVariant.chr) > 0)
				{
					nextEventRegion = false;
				}
				
				// If region is at earlier position on same chromosome, take it
				else if(nextRegion.pos < nextVariant.pos)
				{
					nextEventRegion = true;
				}
				
				// If region is at later position on same chromosome, take variant
				else if(nextRegion.pos > nextVariant.pos)
				{
					nextEventRegion = false;
				}
				
				// Now the case where positions are the same - tie handling depends on overlap mode
				else if(FILTER_MODE.equalsIgnoreCase("CONTAINED_IN_REGION"))
				{
					// Order of priority is variant end, region end, region start, variant start
					if(nextVariant.type == 1)
					{
						nextEventRegion = false;
					}
					else
					{
						nextEventRegion = true;
					}
				}
			}
			
			// After deciding what kind of event to use, process it!
			
			// Case where next event is a region breakpoint
			if(nextEventRegion)
			{
				Event next = regions.get(regionIdx);
				// Start of region
				if(next.type == 1)
				{
					// Add to list of open regions
					openRegions.add(next.id);
				}
				// End of region
				else
				{
					// Remove from list of open regions
					openRegions.remove(next.id);
					
					for(String openVariant : openVariants)
					{
						if(overlaps.get(openVariant).contains(next.id))
						{
							overlaps.get(openVariant).remove(next.id);
						}
					}
				}
				regionIdx++;
			}
			
			// Case where next event is a variant breakpoint
			else
			{
				Event next = variants.get(variantIdx);
				// Start of variant
				if(next.type == 1)
				{
					// Add to list of open variants
					openVariants.add(next.id);
					
					// Initialize list of overlaps to all open regions
					overlaps.put(next.id, new HashSet<String>());
					for(String openRegion : openRegions)
					{
						overlaps.get(next.id).add(openRegion);
					}
				}
				// End of variant
				else
				{
					// Remove from list of open variants
					openVariants.remove(next.id);
					
					// If all overlapping regions were removed, take this out of overlap list
					if(overlaps.get(next.id).size() == 0)
					{
						overlaps.remove(next.id);
					}
				}
				variantIdx++;
			}
		}
		return overlaps;
	}
	
	static void filterVcf() throws Exception
	{
		System.err.println("Getting regions");
		ArrayList<Event> bedEvents = getBedEvents();
		System.err.println("Found " + bedEvents.size() + " region breakpoints");
		System.err.println("Getting variants");
		ArrayList<Event> vcfEvents = getVcfEvents();
		System.err.println("Found " + vcfEvents.size() + " variant breakpoints");
		System.err.println("Finding overlaps");
		HashMap<String, HashSet<String>> overlaps = getOverlaps(bedEvents, vcfEvents);
		System.err.println("Found " + overlaps.size() + " variants with at least one overlap");
		System.err.println("Filtering variants");
		Scanner input = new Scanner(new FileInputStream(new File(vcfFn)));
		PrintWriter out = new PrintWriter(new File(ofn));
		VcfHeader header = new VcfHeader();
		boolean printedHeader = false;
		while(input.hasNext())
		{
			String line = input.nextLine();
			if(line.length() == 0)
			{
				continue;
			}
			if(line.startsWith("#"))
			{
				header.addLine(line);
			}
			else
			{
				if(!printedHeader)
				{
					if(REPORT_MODE.equalsIgnoreCase("INFO"))
					{
						header.addInfoField(reportInfo, "1", "String", "Whether or not the variant is in the regions of interest listed in " + bedFn);
					}
					header.print(out);
					printedHeader = true;
				}
				VcfEntry entry = VcfEntry.fromLine(line);
				String id = entry.getId();
				HashSet<String> curOverlaps = overlaps.getOrDefault(id, null);
				boolean hasOverlap = curOverlaps != null;
				if(entry.getNormalizedType().equals("TRA"))
				{
					String id1 = id + "_breakpoint1", id2 = id + "_breakpoint2";
					HashSet<String> firstOverlap = overlaps.getOrDefault(id1, null);
					HashSet<String> secondOverlap = overlaps.getOrDefault(id2, null);
					hasOverlap = firstOverlap != null && secondOverlap != null;
				}
				if(REPORT_MODE.equalsIgnoreCase("REMOVE"))
				{
					if(hasOverlap)
					{
						out.println(entry);
					}
				}
				if(REPORT_MODE.equalsIgnoreCase("INFO"))
				{
					if(hasOverlap)
					{
						entry.setInfo(reportInfo, "1");
					}
					else
					{
						entry.setInfo(reportInfo, "0");
					}
					out.println(entry);
				}
			}
		}
		input.close();
		out.close();
	}
	
	static class Event implements Comparable<Event>
	{
		String chr;
		int pos;
		int type;
		String id;
		Event(String chr, int pos, int type, String id)
		{
			this.chr = chr;
			this.pos = pos;
			this.type = type;
			this.id = id;
		}
		@Override
		public int compareTo(Event o)
		{
			if(!chr.equals(o.chr))
			{
				return chr.compareTo(o.chr);
			}
			if(pos != o.pos) return pos - o.pos;
			return type - o.type; // Do ends before starts
		}
	}
}
