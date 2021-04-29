import argparse
import bisect
import logging
import sys

import cyvcf2

import utils


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("VCF", type=str)
    parser.add_argument("--supports", default="0,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20")
    parser.add_argument("-o", "--output", default=sys.stdout, type=argparse.FileType("wt"))
    parser.add_argument("--no-out-header", dest="out_header", action="store_false")
    parser.add_argument("--no-out-total-b", dest="out_total_bins", action="store_false")
    parser.add_argument("--no-out-type", dest="out_indiv_types", action="store_false")
    parser.add_argument("--types", type=str, default="INS,DEL,DUP,INV,TRA")
    parser.add_argument("--info-support-field", type=str, default="RE")
    parser.add_argument("--info-reads-field", type=str, default="RNAMES")
    parser.add_argument("--info-type-field", type=str, default="SVTYPE")
    parser.add_argument("--info-len-field", type=str, default="SVLEN")
    args = parser.parse_args()
    logger = logging.getLogger("SV-stats")
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    bins = sorted(set(map(int, args.supports.split(","))))
    logger.debug(f"bins: [{','.join(map(str, bins))}]")
    types = [utils.SVType.from_str(string=s) for s in args.types.split(",")]
    if not set(types).issubset({x for x in utils.SVType}):
        logger.critical(f"Supplied type list {','.join(map(str, types))} is not a subset of a standardized 5 types.")
        exit(1)
    if not all(map(lambda x: x >= 0, bins)):
        logger.warning(f"Some bins were of negative values. Only non-negative values are permitted. Removing all negative values.")
        bins = [x for x in bins if x >= 0]
    if 0 not in bins:
        logger.warning(f"0 value is not in bins. Adding 0.")
        bins = [0] + bins
    bin_counts = {bin_value: {sv_type: 0 for sv_type in types} for bin_value in bins}
    reader = cyvcf2.VCF(args.VCF)
    for cnt, record in enumerate(reader):
        sv_type = utils.get_sv_type(vcf_record=record, info_type_field=args.info_type_field, info_len_field=args.info_len_field, logger=logger)
        sv_support = utils.get_sv_support_cnt(vcf_record=record, info_re_field=args.info_support_field, info_reads_field=args.info_reads_field)
        bin_index = bisect.bisect_right(bins, sv_support)
        assert bin_index > 0
        if sv_type not in bin_counts[bins[bin_index - 1]]:
            continue
        bin_counts[bins[bin_index - 1]][sv_type] += 1
    header = ["bin"]
    if args.out_indiv_types:
        header += types
    if args.out_total_bins:
        header += ["total"]
    if args.out_header:
        print(",".join(map(str, header)), file=args.output)
    type_totals = {sv_type: sum(bin_counts[bin_v][sv_type] for bin_v in bins) for sv_type in types}
    for bin_value in bins:
        sv_type_values = []
        for sv_type in types:
            sv_type_values.append(type_totals[sv_type])
            type_totals[sv_type] -= bin_counts[bin_value][sv_type]
        bin_total = sum(sv_type_values)
        result = f"{bin_value}"
        if args.out_indiv_types:
             result += "," + ",".join(map(str, sv_type_values))
        if args.out_total_bins:
            result += f",{bin_total}"
        print(result, file=args.output)
    assert all(map(lambda x: x == 0, type_totals.values()))


if __name__ == "__main__":
    main()
