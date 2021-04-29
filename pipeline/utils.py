import os
import sys
from collections import defaultdict
from enum import Enum
import logging
from typing import Optional

DEFAULT_THREAD_CNT = 99
DEFAULT_CLUSTER_MEM_MB = 4000

OUTPUT_DIR = "output_dir"
ALIGNMENTS = "alignments"
READS_PATHS = "reads_paths"
SAMPLES = "samples"
LOG = "log"
TOOLS = "tools"
THREADS = "threads"
SAMTOOLS = "samtools"
ALIGNER = "aligner"
MOSDEPTH = "mosdepth"
WINDOW_SIZE = "window_size"
PER_BASE = "per_base"
FAST_MODE = "fast_mode"
PYTHON = "python"

WINNOWMAP = "winnowmap"
MERYL = "meryl"
DISTINCT = "distinct"
K = "k"

AWK = "awk"
SEQTK = "seqtk"
NGMLR = "ngmlr"
PATH = "path"
TMP_DIR = "tmp_dir"
TECH = "tech"
REFERENCE = "ref"
READS_CNT_PER_RUN = "reads_cnt_per_run"


RAW = "raw"
SVS = "svs"
REFINED = "refined"
INS_TO_DUP = "ins_to_dup"
IRIS_REFINED = "iris_refined"
NORM_SV = "norm_sv"
SPECIFIC_MARKED = "specific_marked"
SPEC_READS_FIXED = "spec_reads_fixed"
SPEC_READS_FRACTION = "spec_reads_fraction"
SPEC_LEN = "spec_len"
MAX_DUP_LENGTH = "max_dup_length"

SNIFFLES = "sniffles"
MIN_SUPPORT = "min_support"
MIN_LENGTH = "min_length"
MAX_NUM_SPLIT_READS = "max_num_splits"
MAX_DISTANCE = "max_distance"
NUM_READS_REPORT = "num_reads_report"
MIN_SEQ_SIZE = "min_seq_size"

SV_TOOLS_ENABLED = "sv_tools_enabled"

JAVA = "java"
JASMINE = "jasmine"
IRIS = "iris"
SRC_PATH = "src_path"
MINIMAP2 = "minimap2"
RACON = "racon"
SCRIPT_NAME = "script_name"
MIN_INS_LENGTH = "min_ins_length"
MAX_OUT_LENGTH = "max_out_length"
MAX_INS_DIST = "max_ins_dist"
MAX_LENGTH_CHANGE = "max_length_change"
IS_MERGING = "is_merging"
NORMALIZE_TYPES = "normalize_types"
USE_STRANDS = "use_strands"
USE_TYPES = "use_types"
USE_EDIT_DISTANCE = "use_edit_distance"
USE_END = "use_end"
STRATEGY = "strategy"
KD_TREE_NORM = "kd_tree_norm"
MAX_DISTANCE_LINEAR = "max_distance_linear"
MIN_DISTANCE = "min_distance"
MIN_SEQ_ID = "min_seq_id"
K_JACCARD = "k_jaccard"

SV_SIZES = "sv_sizes"
BINS = "bins"
TYPES = "types"
ABS_LENGTH = "abs_length"
INFO_LENGTH_FIELD = "info_length_field"

MEM_MB_PER_THREAD = "mem_mb_per_thread"
MEM_MB_CORE = "mem_mb_core"

NCPUS = "nCPUs"

STATS = "stats"


ENABLE_SV_INFERENCE = "enable_sv_inference"
ENABLE_SV_REFINEMENT = "enable_sv_refinement"
ENABLE_IS_MERGING = "enable_is_merging"
ENABLE_ALIGNMENT_STATS = "enable_alignment_stats"


EXISTING_ALIGNMENTS = "existing_alignments"
BGZIP = "bgzip"

##########
#
# snakemake data preparation utils
#
##########


def ensure_samples_correctness(config):
    if SAMPLES not in config or not isinstance(config[SAMPLES], list) or len(config[SAMPLES]) < 1:
        raise ValueError("Configuration data file is missing information about samples or the setup is not dictionary-like")


def get_samples_to_reads_paths(config):
    samples_to_reads_paths = defaultdict(list)
    for sample_data in config["samples"]:
        sample_name = sample_data["sample"]
        if TECH not in sample_data or sample_data[TECH].lower() not in ["ont", "pb", "pacbio", "pbccs", "pacbioccs"]:
            raise ValueError(
                f"incorrect or missing tech {sample_data[TECH]} specified for sample {sample_name} in data.yaml. Only ONT or PB are supported, and tech specification is required")
        tech = sample_data[TECH].upper()
        has_alignment = os.path.exists(os.path.join(config.get(OUTPUT_DIR, ""), ALIGNMENTS, f"{sample_name}_{tech}.sort.bam"))
        if not has_alignment:
            if READS_PATHS not in sample_data or not isinstance(sample_data[READS_PATHS], list) or len(sample_data[READS_PATHS]) < 1:
                raise ValueError(
                    f"Error when parsing reads paths for sample {sample_name} sample. Make sure the entries are formatted as a list of strings under the {READS_PATHS} key")
            if (sample_name, tech) in samples_to_reads_paths:
                warning_message = f"sample {sample_name} with read tech {tech} is specified in input data multiple times."
                if not config.get("allow_dup_st_entries", False):
                    raise ValueError(f"Error! {warning_message}")
                else:
                    print(f"WARNING! {warning_message} Proceeding because `allow_dup_st_entries` is set to True", file=sys.stderr)
            for read_path in sample_data[READS_PATHS]:
                if not read_path.endswith(("fastq", "fq", "fastq.gz", "fq.gz", "fasta", "fasta.gz", "fa", "fa.gz")):
                    raise ValueError(f"Unsupported input format for read path {read_path}. Only 'fastq', 'fq', 'fastq.gz', 'fq.gz', 'fasta', 'fasta.gz', 'fa', and 'fa.gz' are supported")
                samples_to_reads_paths[(sample_name, tech)].append(read_path)
            if len(samples_to_reads_paths[(sample_name, tech)]) != len(set(samples_to_reads_paths[(sample_name, tech)])):
                warning_message = f"sample {sample_name} with read tech {tech} has some read file paths specified multiple times."
                if not config.get("allow_dup_reads_entries", False):
                    raise ValueError(f"Error! {warning_message}")
                else:
                    print(f"WARNING! {warning_message} Proceeding because `allow_dup_reads_entries` is set to True", file=sys.stderr)
        else:
            samples_to_reads_paths[(sample_name, tech)].append("")
    return samples_to_reads_paths


def ensure_aligner(config):
    if config['aligner'] not in {"ngmlr", "minimap2", "winnowmap"}:
        raise ValueError(f'unsupported aligner option {config["aligner"]}, only ngmlr, minimap2, and winnowmap are supported')


def get_extra_alignments_paths(config):
    samples_to_reads_paths = defaultdict(list)
    for sample_data in config["samples"]:
        sample_name = sample_data["sample"]
        if TECH not in sample_data or sample_data[TECH].lower() not in ["ont", "pb", "pacbio", "pbccs", "pacbioccs"]:
            raise ValueError(
                f"incorrect or missing tech {sample_data[TECH]} specified for sample {sample_name} in data.yaml. Only ONT or PB are supported, and tech specification is required")
        tech = sample_data[TECH].upper()
        if EXISTING_ALIGNMENTS not in sample_data or not isinstance(sample_data[EXISTING_ALIGNMENTS], list) or len(sample_data[EXISTING_ALIGNMENTS]) < 1:
            samples_to_reads_paths[(sample_name, tech)] = []
            continue
        for alignment_path in sample_data[EXISTING_ALIGNMENTS]:
            if not alignment_path.endswith(("bam")):
                raise ValueError(
                    f"Unsupported extra alignment format for alignment {alignment_path}. Only 'bam' are supported")
            samples_to_reads_paths[(sample_name, tech)].append(alignment_path)
    return samples_to_reads_paths


def get_samples_regex(samples_to_reads_paths):
    return f"({'|'.join(x[0] for x in samples_to_reads_paths.keys())})"


def get_reads_paths_regex(samples_to_reads_paths):
    bases = set()
    for (sample_name, tech), reads_paths in samples_to_reads_paths.items():
        for read_path in reads_paths:
            bases.add(os.path.basename(read_path))
    return f"({'|'.join(bases)})"


def get_tech_regex(config):
    techs = set()
    for sample_data in config[SAMPLES]:
        techs.add(sample_data[TECH])
    return f"({'|'.join(techs)})"


def ensure_ref_correctness(config):
    if REFERENCE not in config:
        raise ValueError(f"No reference fasta file specified under 'ref' key in data.yaml. Reference is required.")


def get_sniffles_sens_suffix(config):
    min_support = config.get(TOOLS, {}).get(SNIFFLES, {}).get(MIN_SUPPORT, 2)
    min_length = config.get(TOOLS, {}).get(SNIFFLES, {}).get(MIN_LENGTH, 20)
    return f"s{min_support}l{min_length}"


SUPPORTED_SV_TOOLS = {"sniffles"}


def ensure_enabled_sv_tools(config):
    for tool in config[SV_TOOLS_ENABLED]:
        if tool.lower() not in SUPPORTED_SV_TOOLS:
            raise ValueError(f"Attempt to enable unsupported SV inference tool {tool}. Only {','.join(SUPPORTED_SV_TOOLS)} are supported")


def get_min_support(coverage_file, min_support_fixed_cnt, min_support_fraction):
    coverage = 100
    with open(coverage_file, "rt") as source:
        for line in source:
            coverage = int(float(line.strip().split("=")[1].strip()))
            print(f"extracted coverage of {coverage} from file {coverage_file}")
            break
    result = min(int(min_support_fixed_cnt), int(coverage * min_support_fraction))
    print(f"target min support cnt {result} with min support fixed cnt = {min_support_fixed_cnt} and min_support_fraction = {min_support_fraction}")
    return result

##########
#
# SV type and length utils
#
##########


class SVType(Enum):
    INS = "INS"
    DEL = "DEL"
    DUP = "DUP"
    INV = "INV"
    TRA = "TRA"

    def __str__(self) -> str:
        return str(self.value)

    def __repr__(self):
        return str(self)

    @classmethod
    def from_str(cls, string: str) -> "SVType":
        for entry in cls:
            if string.lower() == entry.value.lower():
                return entry
        raise ValueError(f"Could not determine SVType from its supplied str version {string}")


def get_chr_from_alt_bnd_record(bnd_string, default: str = "XXX") -> str:
    splitter = "[" if "[" in bnd_string else "]"
    chr_entry = [x for x in bnd_string.split(splitter) if ":" in x]
    if len(chr_entry) < 1:
        return default
    return chr_entry[0].split(":")[0]


def get_sv_type(vcf_record, info_type_field: str = "SVTYPE", info_len_field: str = "SVLEN", logger: Optional[logging.Logger] = None) -> SVType:
    logger = logger if logger else logging.getLogger("Dummy")
    strands = vcf_record.INFO.get("STRANDS", "??")
    chr1 = str(vcf_record.CHROM)
    chr2 = str(vcf_record.INFO.get("CHR2", get_chr_from_alt_bnd_record(bnd_string=vcf_record.ALT[0], default=chr1)))
    if chr1 != chr2:
        return SVType.TRA
    if strands in ["--", "++"]:
        return SVType.INV
    if strands == "-+":
        return SVType.DUP
    info_svtype = vcf_record.INFO.get(info_type_field, None)
    if info_svtype is not None:
        if "INS" in info_svtype:
            return SVType.INS
        if "DEL" in info_svtype:
            return SVType.DEL
    coord_length = get_sv_length_from_coordinates(vcf_record)
    if coord_length in [0, 1]:
        return SVType.INS
    info_length = vcf_record.INFO.get(info_len_field, None)
    if info_length is not None and int(float(info_length)) < 0:
        return SVType.DEL
    logger.warning(f"Can't determine the SV type for VCF record {str(vcf_record)}. Defaulting to DEL")
    return SVType.DEL


def get_sv_length_from_coordinates(vcf_record) -> int:
    try:
        return abs(int(vcf_record.POS) - vcf_record.INFO["END"])
    except KeyError:
        print(f"No END field in VCF record {str(vcf_record)}")


def get_sv_length_from_ref_alt(vcf_record) -> int:
    return abs(len(vcf_record.ALT[0]) - len(vcf_record.REF))


def get_sv_length(vcf_record, abs_value: bool = True, sv_type: Optional[SVType] = None, info_len_field: str = "SVLEN", info_type_field: str = "SVTYPE") -> int:
    """
    0 value is reserved for TRA SVs
    """
    sv_type = sv_type if sv_type else get_sv_type(vcf_record=vcf_record, info_type_field=info_type_field)
    result = 0
    if sv_type == SVType.TRA:
        result = 0
    elif sv_type in [SVType.DUP, SVType.INV]:
        result = int(float(vcf_record.INFO.get(info_len_field, get_sv_length_from_coordinates(vcf_record))))
    elif sv_type == SVType.INS:
        result = int(float(vcf_record.INFO.get(info_len_field, get_sv_length_from_ref_alt(vcf_record))))
    elif sv_type == SVType.DEL:
        result = int(float(vcf_record.INFO.get(info_len_field, get_sv_length_from_coordinates(vcf_record))))
        if result > 0:
            result *= -1
    if abs_value:
        result = abs(result)
    return result


def get_sv_support_cnt(vcf_record, info_re_field: str = "RE", info_reads_field: str = "RNAMES") -> int:
    re_value = int(vcf_record.INFO.get(info_re_field, 0))
    if re_value != 0:
        return re_value
    reads = vcf_record.INFO.get(info_reads_field, "").split(",")
    if len(reads) > 1 or len(reads[0]) > 0:
        return len(reads)
    return 0
