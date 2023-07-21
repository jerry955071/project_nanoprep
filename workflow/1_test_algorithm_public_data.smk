configfile: "configs/0_read_sampling/config.json"

# def query(id:str) -> dict
def query(id:str) -> dict:
    for sample in config["samples"]:
        if sample["id"] == id:
            return sample
    raise ValueError(f"Sample {id} not found in config file")

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
        # perform test1 on fastq samples
        expand(
            "outputs/test_1/{sample_id}",
            sample_id=fastq_sampling_list
        )

# test1
# ==============================================================================
rule test1:
    conda:
        "envs/python-dev_nanoprep_v0.0.13.yaml"
    input:
        "outputs/sampled_fastq/{sample_id}.fq"
    output:
        directory("outputs/test_1/{sample_id}")
    params:
        p5_sense=lambda wildcards: query(wildcards.sample_id)["p5_sense"],
        p3_sense=lambda wildcards: query(wildcards.sample_id)["p3_sense"]
    log:
        "logs/test_1/{sample_id}.log"
    shell:
        """
        python scripts/1_test_algorithm/test1.py \
            {input} \
            {output} \
            {params.p5_sense} \
            {params.p3_sense} \
            > {log} 2>&1
        """
