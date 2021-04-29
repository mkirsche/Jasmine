import os
import utils

if os.path.exists("data.yaml"):
    configfile: "data.yaml"
if os.path.exists("tools.yaml"):
    configfile: "tools.yaml"

output_dir = config.get(utils.OUTPUT_DIR, "")
alignment_output_dir = os.path.join(output_dir, utils.ALIGNMENTS)
svs_output_dir = os.path.join(output_dir, utils.SVS)
raw_svs_output_dir = os.path.join(svs_output_dir, utils.RAW)
refined_svs_output_dir = os.path.join(svs_output_dir, utils.REFINED)
ins_to_dup_output_dir = os.path.join(refined_svs_output_dir, utils.INS_TO_DUP)
normalize_sv_types_output_dir = os.path.join(refined_svs_output_dir, utils.NORM_SV)
iris_refined_output_dir = os.path.join(refined_svs_output_dir, utils.IRIS_REFINED)
specific_marked_output_dir = os.path.join(refined_svs_output_dir, utils.SPECIFIC_MARKED)

sniffles_sens_suffix = utils.get_sniffles_sens_suffix(config)
tech_regex = utils.get_tech_regex(config)
samples_regex = utils.get_samples_regex(sample_to_reads_paths)

samtools_config = config.get(utils.TOOLS, {}).get(utils.SAMTOOLS, {})
java_config=config.get(utils.TOOLS, {}).get(utils.JAVA, {})
jasmine_config=config.get(utils.TOOLS, {}).get(utils.JASMINE, {})
jasmine_is_config=jasmine_config.get(utils.IS_MERGING, {})
iris_config=config.get(utils.TOOLS, {}).get(utils.IRIS, {})
minimap2_config=config.get(utils.TOOLS, {}).get(utils.MINIMAP2, {})
racon_config=config.get(utils.TOOLS, {}).get(utils.RACON, {})
sv_sizes_config=config.get(utils.TOOLS, {}).get(utils.SV_SIZES, {})

rule sv_tally:
    input: os.path.join(refined_svs_output_dir, "{file}")
    output: os.path.join(refined_svs_output_dir, utils.STATS, "{file}.stats.sizes.txt")
    log: os.path.join(refined_svs_output_dir, utils.LOG, "{file}.stats.sizes.txt.log")
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    params:
        python=config.get(utils.TOOLS, {}).get(utils.PYTHON, {}).get(utils.PATH, "python3"),
        script_path=sv_sizes_config.get(utils.PATH, "sv_sizes.py"),
        bins=sv_sizes_config.get(utils.BINS, "1,30,50,100,150,200,350,300,500,750,1000,2000,5000,10000,50000,100000,500000"),
        types=sv_sizes_config.get(utils.TYPES, "INS,DEL,DUP,INV,TRA"),
        abs_length="" if sv_sizes_config.get(utils.ABS_LENGTH, True) else "--no-abs-length",
        info_len_field=sv_sizes_config.get(utils.INFO_LENGTH_FIELD, "SVLEN")
    shell:
         "{params.python} {params.script_path} {input} -o {output} --bins {params.bins} --types {params.types} {params.abs_length} --info-len-field {params.info_len_field} &> {log}"


rule intra_sample_merged_specific:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.{svtypes}ism.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.{svtypes}ism.specific.vcf")
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    shell:
        "awk '($0 ~/^#/ || $0 ~/IS_SPECIFIC=1/)' {input} > {output}"


rule intra_sample_merging:
    input: vcf_list=os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.{svtypes}vcf_list.txt")
    output: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.{svtypes}ism.vcf")
    threads: lambda wc: min(cluster_config.get("intra_sample_merging", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_is_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(refined_svs_output_dir, utils.LOG, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.{svtypes}ism.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 4000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        java_src=jasmine_config.get(utils.SRC_PATH, ""),
        java=java_config.get(utils.PATH, "java"),
        normalize_types="--normalize_type" if jasmine_is_config.get(utils.NORMALIZE_TYPES, True) else "",
        use_types="" if jasmine_is_config.get(utils.USE_TYPES, True) else "--ignore_type",
        use_strands="" if jasmine_is_config.get(utils.USE_STRANDS, True) else "--ignore_strand",
        use_edit_distance="--use_edit_dist" if jasmine_is_config.get(utils.USE_EDIT_DISTANCE, False) else "",
        use_end="--use_end" if jasmine_is_config.get(utils.USE_END, False) else "",
        max_distance=jasmine_is_config.get(utils.MAX_DISTANCE, 200),
        min_distance=jasmine_is_config.get(utils.MIN_DISTANCE, -1),
        strategy="" if jasmine_is_config.get(utils.STRATEGY, "default") == "default" else ("--centroid_merging" if "centroid" in jasmine_is_config.get(utils.STRATEGY, "default") else "--clique_merging"),
        kd_tree_norm=jasmine_is_config.get(utils.KD_TREE_NORM, 2),
        max_dist_linear=jasmine_is_config.get(utils.MAX_DISTANCE_LINEAR, 0),
        min_seq_id=jasmine_is_config.get(utils.MIN_SEQ_ID, 0),
        k_jaccard=jasmine_is_config.get(utils.K_JACCARD, 9)
    shell:  
        "{params.java} -cp {params.java_src} Main file_list={input.vcf_list} {params.normalize_types} {params.use_types} {params.use_strands} "
        "{params.use_edit_distance} {params.strategy} threads={threads} kd_tree_norm={params.kd_tree_norm} max_dist={params.max_distance} --allow_intrasample "
        "max_dist_linear={params.max_dist_linear} min_dist={params.min_distance} min_seq_id={params.min_seq_id} k_jaccard={params.k_jaccard} out_file={output} > {log}"

rule intra_sample_merging_create_vcf_list:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.{svtypes}vcf")
    output: temp(os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.{svtypes}vcf_list.txt"))
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    run:
        dirname = os.path.dirname(output[0])
        os.makedirs(dirname, exist_ok=True)
        with open(output[0], "wt") as dest:
            print(input[0], file=dest)

rule specific_or_sv_types:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.specific.vcf")
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    shell:
        "awk '($0 ~/^#/ || $0 ~/IS_SPECIFIC=1/)' {input} > {output}"

rule spec_marked_sensitive_or_sv_types_ins_to_dup:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.vcf")
    threads: lambda wc: min(cluster_config.get("sensitive_ins_to_dup_conversion", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(refined_svs_output_dir, utils.LOG, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 4000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        ins_to_dup_script_name=jasmine_config.get(utils.SCRIPT_NAME, "InsertionsToDuplications")
    shell:
        "{params.java} -cp {params.java_src} {params.ins_to_dup_script_name} {input} {output} &> {log}"

rule specific_new_sv_types:
    input: os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.specific.vcf")
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    shell:
        "awk '($0 ~/^#/ || $0 ~/IS_SPECIFIC=1/)' {input} > {output}"

rule normalize_sv_types_nSVtypes:
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.norm.vcf")
    input: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.vcf"),
    log: os.path.join(refined_svs_output_dir, utils.LOG, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.norm.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 20000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        tmp_file_path=lambda wc: os.path.join("output", f"{wc.sample}_{wc.tech}_sniffles.{sniffles_sens_suffix}.refined.nSVtypes_normalizeTypes.vcf"),
    shell:
        "{params.java} -cp {params.java_src} Main --preprocess_only --pre_normalize --comma_filelist file_list={input} &> {log} && mv {params.tmp_file_path} {output}"

rule normalize_sv_types_oSVtypes:
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.norm.vcf")
    input: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.vcf"),
    log: os.path.join(refined_svs_output_dir, utils.LOG, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".refined.norm.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 20000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        tmp_file_path=lambda wc: os.path.join("output", f"{wc.sample}_{wc.tech}_sniffles.{sniffles_sens_suffix}.refined_normalizeTypes.vcf"),
    shell:
        "{params.java} -cp {params.java_src} Main --preprocess_only --pre_normalize --comma_filelist file_list={input} &> {log} && mv {params.tmp_file_path} {output}"


rule spec_marked_sensitive_new_sv_type_final_location:
    input: os.path.join(specific_marked_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined_markedSpec.vcf")
    output: os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles." + sniffles_sens_suffix + ".refined.nSVtypes.vcf")
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    params:
        min_support_fixed=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FIXED, 10),
        min_support_fraction=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FRACTION, 0.25),
        min_length=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_LEN, 30),
    run:
        shell("cp {input} {output}")

rule spec_marked_sensitive_new_sv_types:
    input: vcf=os.path.join(iris_refined_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined.vcf"),
           coverage=os.path.join(alignment_output_dir, utils.STATS, "{sample}_{tech}.coverage.txt"),
           vcf_file_list=os.path.join(iris_refined_output_dir, "{sample}_{tech}_sniffles", "vcf_list_dupToIns_irisRefined.txt")
    output: vcf=os.path.join(specific_marked_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined_markedSpec.vcf"),
            vcf_file_list=os.path.join(specific_marked_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list_dupToIns_irisRefined_markedSpec.txt")
    threads: lambda wc: min(cluster_config.get("sensitive_ins_to_dup_conversion", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(specific_marked_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined_markedSpec.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 20000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        output_dir=lambda wc: os.path.join(specific_marked_output_dir, wc.sample + "_" + wc.tech + "_sniffles"),
        min_support_fixed=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FIXED, 10),
        min_support_fraction=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_READS_FRACTION, 0.25),
        min_length=jasmine_config.get(utils.SPECIFIC_MARKED, {}).get(utils.SPEC_LEN, 30),
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
    run:
        min_support=utils.get_min_support(input.coverage, params.min_support_fixed, params.min_support_fraction)
        shell("{params.java} -cp {params.java_src} Main file_list={input.vcf_file_list} --preprocess_only --mark_specific out_dir={params.output_dir} spec_reads=" + str(min_support) + " spec_len={params.min_length} out_file=test.vcf &> {log}")

rule refined_sensitive_new_sv_types:
    input: vcf=os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns.vcf"),
           vcf_file_list=os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "vcf_list_dupToIns.txt"),
           bam=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam"),
           bam_bai=os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam.bai"),
           bam_file_list=os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "bams.txt"),
           ref=config[utils.REFERENCE],
    output: vcf=os.path.join(iris_refined_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined.vcf"),
            vcf_file_list=os.path.join(iris_refined_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list_dupToIns_irisRefined.txt")
    threads: lambda wc: min(cluster_config.get("refined_sensitive_new_sv_types", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), iris_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(iris_refined_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns_irisRefined.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: iris_config.get(utils.MEM_MB_CORE, 20000) + iris_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        output_dir=lambda wc: os.path.join(iris_refined_output_dir, wc.sample + "_" + wc.tech + "_sniffles"),
        iris_output_dir=lambda wc: os.path.join(iris_refined_output_dir, wc.sample + "_" + wc.tech + "_sniffles", "iris"),
        samtools=samtools_config.get(utils.PATH, "samtools"),
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        minimap2=minimap2_config.get(utils.PATH, "minimap2"),
        racon=racon_config.get(utils.PATH, "racon"),
        max_out_length=iris_config.get(utils.MAX_OUT_LENGTH, 100000),
        min_ins_length=iris_config.get(utils.MIN_INS_LENGTH, 30),
        max_ins_dist=iris_config.get(utils.MAX_INS_DIST, 100),
        max_length_change=iris_config.get(utils.MAX_LENGTH_CHANGE, 0.25),
    shell:
        "{params.java} -cp {params.java_src} Main file_list={input.vcf_file_list} --run_iris --preprocess_only genome_file={input.ref} bam_list={input.bam_file_list} "
        "--iris_args=minimap_path={params.minimap2},racon_path={params.racon},samtools_path={params.samtools},threads={threads},out_dir={params.iris_output_dir},"
        "min_ins_length={params.min_ins_length},max_out_length={params.max_out_length},max_ins_dist={params.max_ins_dist},max_length_change={params.max_length_change},--keep_long_variants,"
        "--resume out_dir={params.output_dir} out_file=test.vcf &> {log}"


rule create_bam_file_list:
    output: temp(os.path.join(ins_to_dup_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "bams.txt"))
    input: os.path.join(alignment_output_dir, "{sample}_{tech}.sort.bam")
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    run:
        with open(output[0], "wt") as dest:
            print(input[0], file=dest)


rule sensitive_dup_to_ins_conversion:
    input: vcf=os.path.join(raw_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".vcf"),
           vcf_file_list=os.path.join(refined_svs_output_dir, "{sample}_{tech}_sniffles", "vcf_list.txt"),
           ref=config[utils.REFERENCE],
    output: vcf=temp(os.path.join(ins_to_dup_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns.vcf")),
            vcf_file_list=temp(os.path.join(ins_to_dup_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list_dupToIns.txt"))
    threads: lambda wc: min(cluster_config.get("sensitive_ins_to_dup_conversion", {}).get(utils.NCPUS, utils.DEFAULT_THREAD_CNT), jasmine_config.get(utils.THREADS, utils.DEFAULT_THREAD_CNT))
    log: os.path.join(ins_to_dup_output_dir, "{sample}_{tech}_sniffles", "{sample}_{tech}_sniffles." + sniffles_sens_suffix + "_dupToIns.vcf.log")
    resources:
        mem_mb=lambda wildcards, threads: jasmine_config.get(utils.MEM_MB_CORE, 10000) + jasmine_config.get(utils.MEM_MB_PER_THREAD, 1000) * threads
    params:
        output_dir=lambda wc: os.path.join(ins_to_dup_output_dir, wc.sample + "_" + wc.tech + "_sniffles"),
        java_src=":".join(x for x in [jasmine_config.get(utils.SRC_PATH, ""), iris_config.get(utils.SRC_PATH, "")] if len(x) > 0),
        java=java_config.get(utils.PATH, "java"),
        samtools=samtools_config.get(utils.PATH, "samtools"),
        max_dup_length=jasmine_config.get(utils.INS_TO_DUP, {}).get(utils.MAX_DUP_LENGTH, 10000),
    shell:
         "{params.java} -cp {params.java_src} Main file_list={input.vcf_file_list} --dup_to_ins genome_file={input.ref} "
         "--preprocess_only out_dir={params.output_dir} threads={threads} samtools_path={params.samtools} max_dup_length={params.max_dup_length} out_file=test.vcf &> {log}"

rule create_first_vcf_file_list:
    output: temp(os.path.join(refined_svs_output_dir, "{sample," + samples_regex + "}_{tech," + tech_regex + "}_sniffles", "vcf_list.txt"))
    input: os.path.join(raw_svs_output_dir, "{sample}_{tech}_sniffles." + sniffles_sens_suffix + ".vcf")
    resources:
        mem_mb=utils.DEFAULT_CLUSTER_MEM_MB
    run:
        dirname = os.path.dirname(output[0])
        os.makedirs(dirname, exist_ok=True)
        with open(output[0], "wt") as dest:
            print(input[0], file=dest)

localrules: create_first_vcf_file_list, create_bam_file_list, specific_or_sv_types, specific_new_sv_types, intra_sample_merged_specific, intra_sample_merging_create_vcf_list, sv_tally

include: "call_svs_sniffles_single.snakefile"
