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
        if sample["id"] in id:
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
        fast5_sampling_list.append(sample["id"])

rule all:
    input:
        expand(
            "outputs/transcript_coverage/sampled_fastq/{id}.json",
            id=fastq_sampling_list
        ),
        expand(
            "outputs/transcript_coverage/nanoprep/{id}.{category}.json",
            id=fastq_sampling_list,
            category=["fl", "fu", "tr"]
        ),
        expand(
            "outputs/transcript_coverage/pychopper-{method}/{id}.{category}.json",
            method=["phmm", "edlib"],
            id=fast5_sampling_list,
            category=["fl", "fu", "tr"],
        )
        


# get_transcript_coverage
# ==============================================================================
rule get_transcript_coverage:
    conda: "envs/python-dev_nanoprep_v0.0.15.yaml"
    input: 
        "outputs/minimap2/{name}.sam"
    output:
        "outputs/transcript_coverage/{name}.json"
    log:
        "logs/2_test_nanoprep_public_data/get_transcript_coverage/{name}.log"
    shell:
        """
        ./scripts/2_test_nanoprep/get_transcript_coverage \
            {input} {output} \
            > {log} 2>&1
        """


# minimap2
# ==============================================================================
rule minimap2:
    threads: 10
    input:
        query="outputs/{name}.fq",
        ref=lambda wildcards: query(wildcards.name)["ref_transcripts"]
    output:
        "outputs/minimap2/{name}.sam"
    log:
        "logs/2_test_nanoprep_public_data/minimap2/{name}.log"
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
    threads: 1
    conda: "envs/python-dev_nanoprep_v0.0.15.yaml"
    input:
        "outputs/sampled_fastq/{id}.fq"
    output:
        fl="outputs/nanoprep/{id}.fl.fq",
        fu="outputs/nanoprep/{id}.fu.fq",
        tr="outputs/nanoprep/{id}.tr.fq",
        report="outputs/nanoprep/{id}.report"
    params:
        p5_sense=lambda wildcards: query(wildcards.id)["p5_sense"],
        p3_sense=lambda wildcards: query(wildcards.id)["p3_sense"],
        skip_lowq=7,
        skip_short=0,
        beta=.5,
        n=100000
    log:
        "logs/2_test_nanoprep_public_data/nanoprep/{id}.log"
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
        "outputs/sampled_fastq/{id}.fq"
    output:
        fl="outputs/pychopper-{method}/{id}.fl.fq",
        fu="outputs/pychopper-{method}/{id}.fu.fq",
        tr="outputs/pychopper-{method}/{id}.tr.fq",
        report="outputs/pychopper-{method}/{id}.pdf",
        stats="outputs/pychopper-{method}/{id}.tsv"
    params:
        kit=lambda wildcards: query(wildcards.id)["kit"],
        min_qual=7,
        min_len=0,
        method=lambda wildcards: wildcards.method,
        autotune_nr=100000,
        rescue=lambda wildcards: query(wildcards.id)["kit"]
    log:
        "logs/2_test_nanoprep_public_data/pychopper-{method}/{id}.log"
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