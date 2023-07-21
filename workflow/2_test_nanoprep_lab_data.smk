configfile: "configs/0_read_sampling/config.json"


# configure docker mounting options
docker_mount_opt = ""
for volume in config["volumes"]:
    docker_mount_opt += "-v %s:%s:%s " % (
        volume["real"],
        volume["virtual"],
        volume["mode"]
    )
    if volume["is_workdir"]:
        docker_mount_opt += "-w %s " % volume["virtual"]


# def query(id:str) -> dict
def query(id:str) -> dict:
    for sample in config["samples"]:
        if sample["id"] == id:
            return sample
    raise ValueError(f"Sample {id} not found in config file")


# all
# ==============================================================================
fast5_sampling_list = []
fastq_sampling_list = []
for sample in config["samples"]:
    if "fastq" in sample.keys():
        fastq_sampling_list.append(sample["id"])
    else:
        # skip kit PCS108
        if sample["kit"] == "PCS108":
            continue
        fast5_sampling_list.append(sample["id"])

rule all:
    input:
        expand(
            "outputs/transcript_coverage/sampled_fastq/{guppy_version}/{id}.json",
            guppy_version=config["guppy_versions"],
            id=fast5_sampling_list
        ),
        expand(
            "outputs/transcript_coverage/nanoprep/{guppy_version}/{id}.{category}.json",
            guppy_version=config["guppy_versions"],
            id=fast5_sampling_list,
            category=["fl", "fu", "tr"]
        ),
        expand(
            "outputs/transcript_coverage/pychopper-{method}/{guppy_version}/{id}.{category}.json",
            method=["phmm", "edlib"],
            guppy_version=config["guppy_versions"],
            id=fast5_sampling_list,
            category=["fl", "fu", "tr"],
        )
        

# get_transcript_coverage
# ==============================================================================
rule get_transcript_coverage:
    conda: "envs/python-dev_nanoprep_v0.0.15.yaml"
    input: 
        "outputs/minimap2/{processed}/{guppy_version}/{id}.{category}.sam"
    output:
        "outputs/transcript_coverage/{processed}/{guppy_version}/{id}.{category}.json"
    log:
        "logs/2_test_nanoprep_public_data/get_transcript_coverage/{processed}/{guppy_version}/{id}.{category}.log"
    shell:
        """
        ./scripts/2_test_nanoprep/get_transcript_coverage \
            {input} {output} \
            > {log} 2>&1
        """


# minimap2
# ==============================================================================
rule minimap2:
    threads: 6
    input:
        query="outputs/{processed}/{guppy_version}/{id}.{category}.fq",
        ref=lambda wildcards: query(wildcards.id)["ref_transcripts"]
    output:
        "outputs/minimap2/{processed}/{guppy_version}/{id}.{category}.sam"
    log:
        "logs/2_test_nanoprep_public_data/minimap2/{processed}/{guppy_version}/{id}.{category}.log"
    shell:
        """
        docker run \
            {docker_mount_opt} \
            --rm \
            -u $(id -u) \
            nanozoo/minimap2:2.17--7066fef \
                minimap2 \
                    -t {threads} \
                    -ax map-ont \
                    --eqx \
                    --secondary=no \
                    -o {output} \
                    {input.ref} \
                    {input.query} \
                    > {log} 2>&1
        """


# nanoprep
# ==============================================================================
rule nanoprep:
    threads: 6
    conda: "envs/python-dev_nanoprep_v0.0.17.yaml"
    input:
        "outputs/sampled_fastq/{guppy_version}/{id}.fq"
    output:
        fl="outputs/nanoprep/{guppy_version}/{id}.fl.fq",
        fu="outputs/nanoprep/{guppy_version}/{id}.fu.fq",
        tr="outputs/nanoprep/{guppy_version}/{id}.tr.fq",
        report="outputs/nanoprep/{guppy_version}/{id}.report"
    params:
        p5_sense=lambda wildcards: query(wildcards.id)["p5_sense"],
        p3_sense=lambda wildcards: query(wildcards.id)["p3_sense"],
        skip_lowq=7,
        skip_short=0,
        beta=.5,
        n=100000
    log:
        "logs/2_test_nanoprep_public_data/nanoprep/{guppy_version}/{id}.log"
    shell:
        """
        nanoprep \
            --report {output.report} \
            --beta {params.beta} \
            --p5_sense {params.p5_sense} \
            --p3_sense {params.p3_sense} \
            --skip_lowq {params.skip_lowq} \
            --skip_short {params.skip_short} \
            --input_fq {input} \
            --output_full_length {output.fl} \
            --output_fusion {output.fu} \
            --output_truncated {output.tr} \
            --processes {threads} \
            -n {params.n} \
            > {log} 2>&1
        """


# pychopper
# ==============================================================================
rule pychopper:
    threads: 6
    conda: "envs/nanopore-pychopper_v2.7.6.yaml"
    input:
        "outputs/sampled_fastq/{guppy_version}/{id}.fq"
    output:
        fl="outputs/pychopper-{method}/{guppy_version}/{id}.fl.fq",
        fu="outputs/pychopper-{method}/{guppy_version}/{id}.fu.fq",
        tr="outputs/pychopper-{method}/{guppy_version}/{id}.tr.fq",
        report="outputs/pychopper-{method}/{guppy_version}/{id}.pdf",
        stats="outputs/pychopper-{method}/{guppy_version}/{id}.tsv"
    params:
        kit=lambda wildcards: query(wildcards.id)["kit"],
        min_qual=7,
        min_len=0,
        method=lambda wildcards: wildcards.method,
        autotune_nr=100000,
        rescue=lambda wildcards: query(wildcards.id)["kit"]
    log:
        "logs/2_test_nanoprep_public_data/pychopper-{method}/{guppy_version}/{id}.log"
    shell:
        """
        pychopper \
            -k {params.kit} \
            -Q {params.min_qual} \
            -z {params.min_len} \
            -r {output.report} \
            -u {output.tr} \
            -w {output.fu} \
            -S {output.stats} \
            -Y {params.autotune_nr} \
            -m {params.method} \
            -x {params.rescue} \
            -p \
            -t {threads} \
            {input} \
            {output.fl} \
            > {log} 2>&1
        """

