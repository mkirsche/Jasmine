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
iris_refined_output_dir = os.path.join(refined_svs_output_dir, utils.IRIS_REFINED)
specific_marked_output_dir = os.path.join(refined_svs_output_dir, utils.SPECIFIC_MARKED)

utils.ensure_samples_correctness(config)
sample_to_reads_paths = utils.get_samples_to_reads_paths(config)
utils.ensure_ref_correctness(config)
utils.ensure_enabled_sv_tools(config)


# during development this thing guarantees that only the latest supported part of pipeline produces results
overall_expected_files = []
# print(sample_to_reads_paths)
for (sample, tech) in sample_to_reads_paths.keys():
    if config.get(utils.ENABLE_ALIGNMENT_STATS, True):
        overall_expected_files.append(os.path.join(alignment_output_dir, utils.STATS, f"{sample}_{tech}.coverage.txt"))
        overall_expected_files.append(os.path.join(alignment_output_dir, utils.STATS, f"{sample}_{tech}.samtools.stats.txt"))
        overall_expected_files.append(os.path.join(alignment_output_dir, utils.STATS, f"{sample}_{tech}.alignment.yield.txt"))
        overall_expected_files.append(os.path.join(alignment_output_dir, utils.STATS, f"{sample}_{tech}.mosdepth.global.dist.txt"))
    if config.get(utils.ENABLE_SV_INFERENCE, True):
        for sv_tool in config[utils.SV_TOOLS_ENABLED]:
            if sv_tool == "sniffles":
                suffix = utils.get_sniffles_sens_suffix(config) + "."
            else:
                suffix = ""
            overall_expected_files.append(os.path.join(raw_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}vcf"))
            overall_expected_files.append(os.path.join(raw_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}specific.vcf"))
            overall_expected_files.append(os.path.join(raw_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}vcf.stats.sizes.txt"))
            overall_expected_files.append(os.path.join(raw_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}specific.vcf.stats.sizes.txt"))
            if config.get(utils.ENABLE_SV_REFINEMENT, True):
                overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.specific.vcf"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.specific.vcf.stats.sizes.txt"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.vcf"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.norm.vcf"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.vcf.stats.sizes.txt"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.specific.vcf"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}refined.specific.vcf.stats.sizes.txt"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.vcf"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.norm.vcf"))
                overall_expected_files.append(os.path.join(refined_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}refined.vcf.stats.sizes.txt"))
                if config.get(utils.ENABLE_IS_MERGING, True):
                    overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.ism.vcf"))
                    overall_expected_files.append(os.path.join(refined_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.ism.vcf.stats.sizes.txt"))
                    overall_expected_files.append(os.path.join(refined_svs_output_dir, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.ism.specific.vcf"))
                    overall_expected_files.append(os.path.join(refined_svs_output_dir, utils.STATS, f"{sample}_{tech}_{sv_tool}.{suffix}refined.nSVtypes.ism.specific.vcf.stats.sizes.txt"))

rule main:
    input: overall_expected_files

include: "call_svs_sniffles_single.snakefile"
include: "align_single.snakefile"
include: "jasmine_pre.snakefile"