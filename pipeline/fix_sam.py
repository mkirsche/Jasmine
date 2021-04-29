import argparse
import sys

import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("sam", type=argparse.FileType("rt"), default=sys.stdin)
    parser.add_argument("-o", "--output", type=argparse.FileType("rt"), default=sys.stdout)
    args = parser.parse_args()
    source = pysam.AlignmentFile(args.sam, "r")
    dest = pysam.AlignmentFile(args.output, "w", template=source)
    counter = 0
    while True:
        try:
            record = next(source)
            dest.write(record)
        except OSError as oe:
            print(oe, file=sys.stderr)
            counter += 1
        except StopIteration:
            break
    print(f"{counter} overall skipped records", file=sys.stderr)


if __name__ == "__main__":
    main()
