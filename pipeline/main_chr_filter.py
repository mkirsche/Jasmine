import argparse
import os


def execute(command: str, dry: bool = True):
    if dry:
        print(command)
    else:
        os.system(command)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("VCF")
    parser.add_argument("--chr_bed", default="/work-zfs/mschatz1/pipelines/resources/ref/human/GRCh38.main_chr.bed")
    parser.add_argument("--bad-chr-FLs", default="LKG")
    parser.add_argument("-n", "--dry-run", dest="dry", action="store_true")
    args = parser.parse_args()
    basename = os.path.basename(args.VCF)
    new_vcf_name = basename.split(".")[0] + ".all_chr." + ".".join(basename.split(".")[1:])
    execute(f'mv {args.VCF} {new_vcf_name}', dry=args.dry)
    execute(f'grep "#" {new_vcf_name} > {args.VCF}', dry=args.dry)
    execute(f'bedtools intersect -a {new_vcf_name} -b {args.chr_bed} -u | awk \'{{if ($8 !~/SVTYPE=BND/ || $5 !~ /[\\[\\]][{args.bad_chr_FLs}]/) print $0}}\' >> {args.VCF}',
            dry=args.dry)


if __name__ == "__main__":
    main()
