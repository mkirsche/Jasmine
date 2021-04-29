import argparse
import bisect
import logging
import sys
import cyvcf2

import utils


def bin_string_repr(bin_value: int) -> str:
    suffixes = ["", "K", "M", "G"]
    suffix_index = 0
    while abs(bin_value) >= 1000:
        suffix_index += 1
        bin_value /= 1000
    if bin_value == int(bin_value):
        bin_value = int(bin_value)
    return f"{bin_value:,}{suffixes[suffix_index]}"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("VCF", type=str)
    parser.add_argument("--bins", default="1,30,50,100,150,200,350,300,500,750,1000,2000,5000,10000,50000,100000,500000")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("--no-out-header", dest="out_header", action="store_false")
    parser.add_argument("--no-out-total-t", dest="out_total_type", action="store_false")
    parser.add_argument("--no-out-total-s", dest="out_total_size", action="store_false")
    parser.add_argument("--no-out-indiv-bins", dest="out_indiv_bins", action="store_false")
    parser.add_argument("--types", type=str, default="INS,DEL,DUP,INV,TRA")
    parser.add_argument("--no-abs-length", dest="abs_length", action="store_false")
    parser.add_argument("--info-len-field", type=str, default="SVLEN")
    parser.add_argument("--info-type-field", type=str, default="SVTYPE")
    args = parser.parse_args()
    logger = logging.getLogger("SV-stats")
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    supplied_bins = sorted(set(map(int, args.bins.split(","))))
    bins = [-3000000000] + supplied_bins + [3000000000]
    logger.debug(f"bins: [{','.join(map(str, bins))}]")
    types = [utils.SVType.from_str(string=s) for s in args.types.split(",")]
    if not set(types).issubset({x for x in utils.SVType}):
        logger.critical(f"Supplied type list {','.join(map(str, types))} is not a subset of a standardized 5 types.")
        exit(1)
    bin_counts = {bin_value: {sv_type: 0 for sv_type in types} for bin_value in bins}
    reader = cyvcf2.VCF(args.VCF)
    for cnt, record in enumerate(reader):
        sv_type = utils.get_sv_type(vcf_record=record, info_type_field=args.info_type_field, info_len_field=args.info_len_field, logger=logger)
        sv_length = utils.get_sv_length(record, sv_type=sv_type, abs_value=args.abs_length, info_len_field=args.info_len_field, info_type_field=args.info_type_field)
        bin_index = bisect.bisect_right(bins, sv_length)
        if bin_index < 1:
            logger.error(f"Something is wrong with length bin determination for record {str(record)} with type {str(sv_type)}")
        if sv_type not in bin_counts[bins[bin_index - 1]]:
            continue
        bin_counts[bins[bin_index - 1]][sv_type] += 1
    type_totals = {sv_type: sum(bin_counts[bin_v][sv_type] for bin_v in bins) for sv_type in types}
    header = ["bin"] + types
    if args.out_total_size:
        header += ["total"]
    if args.out_header:
        print(",".join(map(str, header)), file=args.output)
    if args.out_indiv_bins:
        for lv, rv in zip(bins[:-1], bins[1:]):
            bin_str_value = f"[{bin_string_repr(lv)} - {bin_string_repr(rv)})"
            sv_type_values = [bin_counts[lv][sv_type] for sv_type in types]
            bin_total = sum(sv_type_values)
            result = f"{bin_str_value}," + ",".join(map(str, sv_type_values))
            if args.out_total_size:
                result += f",{bin_total}"
            print(result, file=args.output)
    if args.out_total_type:
        print("total," + ",".join(map(str, (type_totals[sv_type] for sv_type in types))), file=args.output)


if __name__ == "__main__":
    main()
