import os
import utils
from collections import defaultdict

if os.path.exists("data.yaml"):
    configfile: "data.yaml"
if os.path.exists("tools.yaml"):
    configfile: "tools.yaml"


output_dir = config.get(utils.OUTPUT_DIR, "")
alignment_output_dir = os.path.join(output_dir, utils.ALIGNMENTS)

utils.ensure_samples_correctness(config)
samples_to_reads_paths = utils.get_samples_to_reads_paths(config)
samples_to_extra_alignments_paths = utils.get_extra_alignments_paths(config)
samples_to_basename_readpaths = defaultdict(dict)
for (sample_name, tech), read_paths in samples_to_reads_paths.items():
    for read_path in read_paths:
        samples_to_basename_readpaths[(sample_name, tech.upper())][os.path.basename(read_path)] = read_path
utils.ensure_ref_correctness(config)

ngmlr_config = config.get(utils.TOOLS, {}).get(utils.NGMLR, {})
samtools_config = config.get(utils.TOOLS, {}).get(utils.SAMTOOLS, {})
mosdepth_config = config.get(utils.TOOLS, {}).get(utils.MOSDEPTH, {})
sed_config = config.get(utils.TOOLS, {}).get("sed", {})
awk_config = config.get(utils.TOOLS, {}).get(utils.AWK, {})
minimap2_config=config.get(utils.TOOLS, {}).get(utils.MINIMAP2, {})
seqtk_config = config.get(utils.TOOLS, {}).get(utils.SEQTK, {})
winnowmap_config=config.get(utils.TOOLS, {}).get(utils.WINNOWMAP, {})
meryl_config=config.get(utils.TOOLS, {}).get(utils.MERYL, {})

samples_regex = utils.get_samples_regex(samples_to_reads_paths)
read_paths_regex = utils.get_reads_paths_regex(samples_to_reads_paths)
tech_regex = utils.get_tech_regex(config)

def split_fastx_dirs(wildcards):
    extensions = read_extensions_per_sample(sample=wildcards.sample, tech=wildcards.tech)
    result = []
    if "fasta" in extensions:
        result.append(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fasta"))
    if "fastq" in extensions:
        result.append(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fastq"))
    return result

rule merged_coverage_mosdepth:
    input: bam=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam"),
           index=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam.bai"),
    output: os.path.join(alignment_output_dir, utils.STATS, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.mosdepth.global.dist.txt")
    message: "Computing mosdepth coverage stats on {input}"
    threads: lambda wildcards: min(cluster_config.get("merged_coverage_mosdepth", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), mosdepth_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.mosdepth.txt.log")
    resources:
        mem_mb=lambda wildcards, threads: mosdepth_config.get(utils.MEM_MB_CORE, 2000) + mosdepth_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        mosdepth=mosdepth_config.get(utils.PATH, "mosdepth"),
        per_base=lambda wc: "" if bool(mosdepth_config.get(utils.PER_BASE, False)) else "-n",
        fast_mode=lambda wc: "--fast-mode" if bool(mosdepth_config.get(utils.FAST_MODE, True)) else "",
        window_size=mosdepth_config.get(utils.WINDOW_SIZE, 500),
        prefix=lambda wc: os.path.join(alignment_output_dir, utils.STATS, f"{wc.sample}_{wc.tech}"),
    shell:
        "{params.mosdepth} {params.per_base} {params.fast_mode} --by {params.window_size} -t {threads} {params.prefix} {input.bam}"


rule merged_average_coverage_samtools:
    input: bam=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam"),
    output: os.path.join(alignment_output_dir, utils.STATS, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.coverage.txt")
    message: "Computing average alignment read depth coverage on {input}"
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.coverage.txt.log")
    resources:
        mem_mb=lambda wildcards, threads: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        samtools=samtools_config.get(utils.PATH, "samtools"),
        awk=awk_config.get(utils.PATH, "awk")
    shell:
        "{params.samtools} depth -a {input.bam} | {params.awk} \'{{sum += $3}} END {{print \"Average coverage (all) = \",sum/NR}}\' > {output} 2> {log}"

rule alignment_yield:
    input: bam=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam")
    output: os.path.join(alignment_output_dir, utils.STATS, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.alignment.yield.txt")
    message: "Computing total read length (i.e., sequencing yield) from the produced bam file."
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.alignment.yield.log")
    resources:
        mem_mb=lambda wildcards, threads: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        samtools=samtools_config.get(utils.PATH, "samtools"),
        awk=awk_config.get(utils.PATH, "awk")
    shell:
        "{params.samtools} bam2fq {input} | awk \'NR%4==2 {{l+=length($0)}} END {{print \"Total read length \",l}}\' > {output} 2> {log}"

rule samtools_stats:
    output: os.path.join(alignment_output_dir, utils.STATS, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.samtools.stats.txt")
    input: bam=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam"),
           ref=config[utils.REFERENCE],
    message: "Computing samtools stats from the produced bam file"
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.samtools.stats.log")
    threads: lambda wildcards: min(cluster_config.get("single_sam_to_sort_bam", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), samtools_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    resources:
        mem_mb=lambda wildcards, threads: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        samtools=samtools_config.get(utils.PATH, "samtools"),
    shell:
        "{params.samtools} stats -r {input.ref} {input.bam} > {output} 2> {log}"


def read_extensions_per_sample(sample, tech):
    result = set()
    for read_path in samples_to_reads_paths[(sample, tech.upper())]:
        if read_path.endswith(("fastq", "fq", "fastq.gz", "fq.gz")):
            result.add("fastq")
        elif read_path.endswith(("fasta", "fa", "fasta.gz", "fa.gz")):
            result.add("fasta")
    return sorted(result)


def aggregated_input_for_bam_merging(wildcards):
    extensions = read_extensions_per_sample(sample=wildcards.sample, tech=wildcards.tech)
    # assert "fasta" in extensions or "fastq" in extensions
    result = []
    if "fasta" in extensions:
        chekpoint_output = checkpoints.split_fasta.get(**wildcards).output[0]
        result.extend(expand(
            os.path.join(alignment_output_dir, f"{wildcards.sample}_{wildcards.tech}_fasta_" + "{chunk_id}.sort.bam"),
            chunk_id=glob_wildcards(os.path.join(chekpoint_output, f"{wildcards.sample}_{wildcards.tech}_fasta_" + "{chunk_id}.fasta.gz")).chunk_id
        ))
    if "fastq" in extensions:
        chekpoint_output = checkpoints.split_fastq.get(**wildcards).output[0]
        result.extend(expand(
            os.path.join(alignment_output_dir, f"{wildcards.sample}_{wildcards.tech}_fastq_" + "{chunk_id}.sort.bam"),
            chunk_id=glob_wildcards(os.path.join(chekpoint_output, f"{wildcards.sample}_{wildcards.tech}_fastq_" + "{chunk_id}.fastq.gz")).chunk_id
        ))
    return result

rule single_bam_to_sort_bam:
    output:temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}_{chunk_id,[a-z]+}.sort.bam"))
    input:
        bam=os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}_{chunk_id}.bam"),
        reads=os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}.{seq_format}.gz"),
        tmp_dir=lambda wc: os.path.join(config["tools"].get(utils.TMP_DIR, ""), f"samtools_tmp_{wc.sample}_{wc.tech}_{wc.seq_format}_{wc.chunk_id}")
    threads:lambda wildcards: min(cluster_config.get("single_sam_to_sort_bam", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), samtools_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Transforming an alignment sam file {input.bam} into a sorted bam file {output}. Requested mem {resources.mem_mb}M on {threads} threads. Cluster config "
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}.sort.bam.log")
    resources:
        mem_mb=lambda wildcards, threads: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        tmp_dir=lambda wc: os.path.join(config["tools"].get(utils.TMP_DIR, ""), f"samtools_tmp_{wc.sample}_{wc.tech}_{wc.seq_format}_{wc.chunk_id}" + os.path.sep),
        samtools=samtools_config.get(utils.PATH, "samtools"),
        mem_mb_per_thread=samtools_config.get(utils.MEM_MB_PER_THREAD, 1000),
    shell:
        "{params.samtools} sort -O bam -o {output} -@ {threads} -m {params.mem_mb_per_thread}M -T {params.tmp_dir} {input.bam} &> {log}"

rule clear_corrupt_sam_alignments:
    output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}" + "_{chunk_id,[a-z]+}.fixed.sam"))
    input: os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format,(fastq|fasta)}" + "_{chunk_id}.sam")
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}" + "_{chunk_id,[a-z]+}.fixed.sam.log")
    threads: 1
    resources:
        mem_mb=lambda wildcards, threads: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        command=lambda wc: f"python {config.get('tools', {}).get('sam_fix', {}).get('path', 'fix_sam.py')}" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "ngmlr" else "cat"
    shell:
        "{params.command} {input} > {output} 2> {log}"


def get_aligner(sample=None, tech=None, default="ngmlr"):
    if sample is None or tech is None:
        return default
    result = default
    for sample_data in config["samples"]:
        if sample_data["sample"] == sample and sample_data["tech"] == tech:
            result = sample_data.get(utils.ALIGNER, result)
            break
    return config.get(utils.ALIGNER, result).lower()



aligners_configs = {
    "ngmlr": ngmlr_config,
    "minimap2": minimap2_config,
    "winnowmap": winnowmap_config,
}

def get_aligner_path(aligner, sample=None, tech=None, default=None):
    if default is None:
        default = aligner
    if sample is None or tech is None:
        return default
    result = default
    result = aligners_configs[aligner].get(utils.PATH, result)
    result = aligners_configs[aligner].get(tech, {}).get(utils.PATH, result)
    for sample_data in config["samples"]:
        if sample_data["sample"] == sample and sample_data["tech"] == tech:
            result = sample_data.get(aligner, {}).get(utils.PATH, result)
            break
    return result

def get_aligner_preset(aligner, tech):
    if "ont" in tech.lower():
        return "map-ont"
    else:
        if aligner == "ngmlr":
            return "map-pacbio"
        elif aligner == "minimap2":
            if "ccs" in tech.lower():
                return "asm20"
            else:
                return "map-pb"
        else:
            if "ccs" in tech.lower():
                return "map-pb"
            else:
                return "map-pb-clr"



rule single_alignment:
    output:temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}" + "_{chunk_id,[a-z]+}.bam"))
    input: reads=os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}.{seq_format}.gz"),
           meryld_db=lambda wc: f"{config[utils.REFERENCE]}_k{meryl_config.get(utils.K, 15)}.txt" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "winnowmap" else os.path.join(alignment_output_dir, f"{wc.sample}_{wc.tech}_{wc.seq_format}", f"{wc.sample}_{wc.tech}_{wc.seq_format}_{wc.chunk_id}.{wc.seq_format}.gz"),
           reference=config[utils.REFERENCE],
    threads: lambda wildcards: min(cluster_config.get("single_alignment", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), aligners_configs[get_aligner(sample=wildcards.sample, tech=wildcards.tech)].get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    message: "Aligning reads from {input} to {output}. Requested mem {resources.mem_mb}M on {threads} threads. Cluster config "
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}.sam.log")
    resources:
        mem_mb=lambda wildcards, threads: aligners_configs[get_aligner(sample=wildcards.sample, tech=wildcards.tech)].get(utils.MEM_MB_CORE, 6000) + aligners_configs[get_aligner(sample=wildcards.sample, tech=wildcards.tech)].get(utils.MEM_MB_PER_THREAD, 1000) * threads,
    params:
        aligner=lambda wc: get_aligner_path(aligner=get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr"), sample=wc.sample, tech=wc.tech),
        input_flag=lambda wc: "-q" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "ngmlr" else "",
        preset_value=lambda wc: get_aligner_preset(aligner=get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr"), tech=wc.tech),
        sam_output_flag=lambda wc: "" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "ngmlr" else "-a",
        reference_flag=lambda wc: "-r" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "ngmlr" else "",
        bam_fix_flag=lambda wc: "--bam-fix" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "ngmlr" else "",
        md_flag=lambda wc: "" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "ngmlr" else "--MD",
        w_db_flag=lambda wc: f"-W {config[utils.REFERENCE]}_k{meryl_config.get(utils.K, 15)}" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "winnowmap" else "",
        w_flag=lambda wc: f"-w {winnowmap_config.get('w', 50)}" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "winnowmap" else "",
        k_flag=lambda wc: f"-k {meryl_config.get(utils.K, 15)}" if get_aligner(sample=wc.sample, tech=wc.tech, default="ngmlr") == "winnowmap" else "",
        samtools=samtools_config.get(utils.PATH,"samtools"),
    shell:
         "{params.aligner} {params.w_db_flag} {params.w_flag} {params.k_flag} {params.reference_flag} {input.reference} {params.input_flag} {input.reads} -t {threads}"
         " -x {params.preset_value} {params.sam_output_flag} {params.bam_fix_flag} {params.md_flag} 2> {log} | {params.samtools} view -O bam -o {output} -"

rule meryl_db_repetitive_extract:
    output: config[utils.REFERENCE] + "_k{k,\d+}.txt"
    input: config[utils.REFERENCE] + "_k{k}DB"
    log: config[utils.REFERENCE] + "_k{k}.log"
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    params:
        meryl=meryl_config.get(utils.PATH, "meryl"),
        distinct=meryl_config.get(utils.DISTINCT, 0.9998),
    shell:
        "{params.meryl} print greater-than distinct={params.distinct} {input} > {output} 2> {log}"


rule meryl_db_creation:
    output: directory(config[utils.REFERENCE] + "_k{k,\d+}DB")
    input: config[utils.REFERENCE]
    log: config[utils.REFERENCE] + "_k{k}DB.log"
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    params:
        meryl=meryl_config.get(utils.PATH, "meryl"),
        k=lambda wc: wc.k,
    shell:
        "{params.meryl} count k={params.k} output {output} {input} &> {log}"


# rule ensure_reads_input_extension:
#     input: os.path.join(alignment_output_dir, "{sample}_{tech}_{seq_format}", "{sample}_{tech}_{seq_format}_{chunk_id}")
#     output: temp(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_{seq_format,(fastq|fasta)}_{chunk_id,[a-z]+}.{seq_format}"))
#     shell: "mv {input} {output} && touch {input}"


def get_fastx_files(sample, tech, extension):
    result = []
    for read_path in samples_to_reads_paths[(sample, tech.upper())]:
        if read_path.endswith(extension):
            result.append(read_path)
    return result


def get_reads_batch_cnt(aligner, default=200000):
    a_config = aligners_configs[aligner]
    return a_config.get(utils.READS_CNT_PER_RUN, default)

checkpoint split_fastq:
    output:
        temp(directory(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fastq")))
    input:
        fastq=lambda wc: get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fastq", "fq")),
        fastq_gz=lambda wc: get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fastq.gz", "fa.gz")),
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    params:
        cat_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fastq", "fq"))) == 0 else f"<(cat {' '.join(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=('fastq', 'fq')))})",
        zcat_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fastq.gz", "fq.gz"))) == 0 else f"<(zcat {' '.join(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=('fastq.gz', 'fq.gz')))})",
        prefix=lambda wc: os.path.join(alignment_output_dir, f"{wc.sample}_{wc.tech}_fastq", f"{wc.sample}_{wc.tech}_fastq_"),
        fastq_cnt=lambda wc: get_reads_batch_cnt(aligner=get_aligner(sample=wc.sample, tech=wc.tech), default=200000) * 4,
        bgzip=config.get(utils.TOOLS,{}).get(utils.BGZIP,{}).get(utils.PATH,"bgzip"),
    shell:
        "mkdir -p {output} && cat {params.cat_command} {params.zcat_command} |"
        " split -l {params.fastq_cnt} -a 3 --filter='{params.bgzip} -c > $FILE.fastq.gz' - {params.prefix}"

checkpoint split_fasta:
    output:
        temp(directory(os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_fasta")))
    input:
        fasta=lambda wc: get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fasta", "fa")),
        fasta_gz=lambda wc: get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fasta.gz", "fa.gz")),
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    params:
        cat_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fasta", "fa"))) == 0 else "<(cat {fasta})".format(fasta=" ".join(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fasta", "fa")))),
        zcat_command=lambda wc: "" if len(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fasta.gz", "fa.gz"))) == 0 else "<(zcat {fasta_gz})".format(fasta_gz=" ".join(get_fastx_files(sample=wc.sample, tech=wc.tech, extension=("fasta.gz", "fa.gz")))),
        prefix=lambda wc: os.path.join(alignment_output_dir, f"{wc.sample}_{wc.tech}_fasta", f"{wc.sample}_{wc.tech}_fasta_"),
        fasta_cnt=lambda wc: get_reads_batch_cnt(aligner=get_aligner(sample=wc.sample, tech=wc.tech), default=200000) * 2,
        seqtk_command=lambda wc: seqtk_config.get(utils.PATH, "seqtk"),
        bgzip=config.get(utils.TOOLS, {}).get(utils.BGZIP, {}).get(utils.PATH, "bgzip"),
    shell:
        "mkdir -p {output} && cat {params.cat_command} {params.zcat_command} |"
        " seqtk seq - |"
        " split -l {params.fasta_cnt} -a 3 --filter='{params.bgzip} -c > $FILE.fasta.gz' - {params.prefix}"

rule samtools_tmp_dir:
    output: temp(directory(os.path.join(config["tools"].get(utils.TMP_DIR, ""), "samtools_tmp_{sample}_{tech}_{seq_format}_{chunk_id}")))
    shell: "mkdir -p {output}"


def get_aggregate_dirs(sample, tech):
    extensions = read_extensions_per_sample(sample=sample, tech=tech)
    return [os.path.join(alignment_output_dir, f"{sample}_{tech}_{ex}") for ex in ["fasta", "fastq"] if ex in extensions]

rule index_bam:
    input: os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam")
    output: os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.sort.bam.bai")
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.sort.bam.bai.log")
    resources:
        mem_mb=lambda wildcards: samtools_config.get(utils.MEM_MB_CORE, 2000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000)
    params:
        samtools=samtools_config.get(utils.PATH, "samtools"),
    shell:
        "{params.samtools} index {input}"

rule merge_sorted:
    output: os.path.join(alignment_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}.sort.bam")
    input: created_bams=aggregated_input_for_bam_merging,
           existing_alignments=lambda wc: samples_to_extra_alignments_paths[(wc.sample, wc.tech)],
           agg_dirs=lambda wc: get_aggregate_dirs(sample=wc.sample, tech=wc.tech)
    message: "Combining sorted bam files. Requested mem {resources.mem_mb}M."
    log: os.path.join(alignment_output_dir, utils.LOG, "{sample}_{tech}.sort.bam.log")
    resources:
        mem_mb=lambda wildcards: samtools_config.get(utils.MEM_MB_CORE, 5000) + samtools_config.get(utils.MEM_MB_PER_THREAD, 1000)
    params:
        samtools=samtools_config.get(utils.PATH, "samtools"),
    shell:
        "{params.samtools} merge {output} {input.created_bams} {input.existing_alignments} &> {log}"


localrules: ensure_reads_input_extension, samtools_tmp_dir
