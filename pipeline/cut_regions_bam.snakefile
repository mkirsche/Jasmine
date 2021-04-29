import os

files = []
out_dir = config["out_dir"]
suffix = config.get("name", "regions")
input_by_base = {}
for path in config["bams"]:
    basename = os.path.basename(path)
    base = os.path.splitext(basename)[0]
    input_by_base[base] = path
    files.append(os.path.join(out_dir, f"{base}.{suffix}.sort.bam"))
    files.append(os.path.join(out_dir, f"{base}.{suffix}.sort.bam.bai"))

regions = []
with open(config["regions"], "rt") as source:
    for line in source:
        data = line.strip().split("\t")
        regions.append(f"{data[0]}:{data[1]}-{data[2]}")
regions = " ".join(regions)


rule all:
    input: files


rule index_bam:
    output: os.path.join(out_dir, "{base}.{suffix," + suffix + "}.sort.bam.bai")
    input: os.path.join(out_dir, "{base}.{suffix}.sort.bam")
    shell:
        "samtools index {input}"

rule sort_cut_bam:
    output: os.path.join(out_dir, "{base}.{suffix," + suffix+ "}.sort.bam")
    input: os.path.join(out_dir, "{base}.{suffix}.bam")
    shell:
        "samtools sort -@ 4 -O bam -o {output} {input}"

rule create_cut_bam:
    output: temp(os.path.join(out_dir, "{base}.{suffix," + suffix + "}.bam"))
    input: bam=lambda wc: input_by_base[wc.base]
    params:
        regions=regions,
    shell:
        "samtools view -O bam -o {output} {input} {params.regions}"