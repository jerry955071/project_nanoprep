configfile: "configs/0_read_sampling/config.json"

# set random seed
RANDOM_SEED = 42

# guppy_basecall_server address
guppy_server = "guppy_basecall_server-basecall_server-1:5566"

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
fast5_sampling_list = []
fastq_sampling_list = []
for sample in config["samples"]:
    if "fastq" in sample.keys():
        fastq_sampling_list.append(sample["id"])
    else:
        fast5_sampling_list.append(sample["id"])


rule all:
    input:
        # perform basecalling for fast5 samples
        expand(
            "outputs/sampled_fastq/{guppy_version}/{sample_id}.fq",
            guppy_version=config["guppy_versions"],
            sample_id=fast5_sampling_list
        ),
        # perform sampling for fastq samples
        expand(
            "outputs/sampled_fastq/{sample_id}.fq",
            sample_id=fastq_sampling_list
        )


# aggregate basecalling results
# ==============================================================================
rule aggregate:
    input:
        "outputs/basecalling/{guppy_version}/{sample_id}"
    output:
        "outputs/sampled_fastq/{guppy_version}/{sample_id}.fq"
    log:  
        "logs/0_read_sampling/aggregate/{guppy_version}/{sample_id}.log"
    shell:
        """
        echo $(date +%Y%m%d-%H:%M:%S) > {log}
        for fastq in $(find {input} -name "*.fastq"); do
            cat $fastq >> {output}
        done
        """


# basecalling
# ==============================================================================
rule basecalling:
    resources:
        gpu=1
    input:
        "outputs/sampled_fast5/{sample_id}/"
    output:
        save_path=directory("outputs/basecalling/{guppy_version}/{sample_id}")
    log:
        "logs/0_read_sampling/basecalling/{guppy_version}/{sample_id}.log"
    params:
        remote_host="b05b01002@192.168.50.79",
        remote_docker="docker-olaf",
        guppy_image=\
            lambda wildcards: 
                f"genomicpariscentre/guppy-gpu:{wildcards.guppy_version}",
        config=\
            lambda wildcards: 
                query(wildcards.sample_id)["config"],
        tmpdir="/home/b05b01002/HDD/guppy_basecaller_tmp/"
    shell:
        """
        ./scripts/basecaller_remote \
            {params.remote_host} \
            {params.remote_docker} \
            {params.guppy_image} \
            {input} \
            {output.save_path} \
            {params.config} \
            {params.tmpdir} \
            > {log} 2>&1
        """


# fast5 sampling
# ==============================================================================
rule fast5_sampling:
    threads: 8
    input:
        fast5_dir = lambda wildcards: query(id=wildcards["sample_id"])["fast5_dir"],
        read_id_list = lambda wildcards: query(id=wildcards["sample_id"])["read_id_list"]
    output:
        sampled_names = temp("outputs/sampled_fast5/{sample_id}.txt"),
        out_dir = directory("outputs/sampled_fast5/{sample_id}/")
    log:
        "logs/0_read_sampling/fast5_sampling/{sample_id}.log"
    shell:
        """
        # echo start time to log file
        echo $(date +%Y%m%d-%H:%M:%S) > {log}
        
        # check number of columns in read_id_list
        NF=$(cat {input.read_id_list} | awk 'NR == 1 {{print NF}}')
        if [[ $NF != 1 ]]; then
            # write header for multi-column read_id_list
            head -n 1 {input.read_id_list} > {output.sampled_names}

            # remove header line
            tail -n +2 {input.read_id_list} > .tmp.read_id_list.{wildcards.sample_id}.txt
        
            # sample 100k reads (append to header)
            ./scripts/random-lines {input.read_id_list} - 100000 {RANDOM_SEED} >> {output.sampled_names}
                
            # clean up tmp file
            rm .tmp.read_id_list.{wildcards.sample_id}.txt
        else
            # sample 100k reads
            ./scripts/random-lines {input.read_id_list} {output.sampled_names} 100000 {RANDOM_SEED}
        fi
        
        # create output directory
        mkdir -p {output.out_dir}

        # extract fast5 records using `fast5_subset`
        docker run \
            {docker_mount_opt} \
            -u $(id -u) \
            --rm \
            --name {wildcards.sample_id}_fast5_subset \
            ccc/ont-fast5-api:latest \
            fast5_subset \
                -i {input.fast5_dir} \
                -s {output.out_dir} \
                -l {output.sampled_names} \
                -n 100000 \
                -r \
                -t {threads} \
            >> {log} 2>&1
        """


# fastq sampling
# ==============================================================================
rule fastq_sampling:
    input:
        lambda wildcards: query(id=wildcards["sample_id"])["fastq"]
    output:
        "outputs/sampled_fastq/{sample_id}.fq"
    log:
        "logs/0_read_sampling/fastq_sampling/{sample_id}.log"
    shell:
        """
        # echo start time to log file
        echo $(date +%Y%m%d-%H:%M:%S) > {log}

        # sample 100k reads
        ./scripts/fastq-subset {input} {output} 100000 {RANDOM_SEED} >> {log} 2>&1 
        """


# UNUSED RULES
# ==============================================================================
# rule sampling_ptr_108_bio1:
#     threads: 10
#     input:
#         filelist = lambda wildcards: query(id="ptr_108_bio1")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/ptr_108_bio1.txt"),
#         out_dir = directory("outputs/sampled_fast5/ptr_108_bio1/")
#     log:
#         "logs/0_read_sampling/sampling_ptr_108_bio1.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         shuf -n 100000 {input.filelist} > {output.sampled_names}
#         mkdir -p {output.out_dir}
#         fast5_dir=$(dirname {input.filelist})
#         cat {output.sampled_names} \
#             | xargs -I % -P {threads} cp $fast5_dir/% {output.out_dir}
#         """
    
# rule sampling_ptr_108_bio2:
#     threads: 10
#     input: 
#         filelist = lambda wildcards: query(id="ptr_108_bio2")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/ptr_108_bio2.txt"),
#         out_dir = directory("outputs/sampled_fast5/ptr_108_bio2/")
#     log:
#         "logs/0_read_sampling/sampling_ptr_108_bio2.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         shuf -n 100000 {input.filelist} > {output.sampled_names}
#         mkdir -p {output.out_dir}
#         fast5_dir=$(dirname {input.filelist})
#         cat {output.sampled_names} \
#             | xargs -I % -P {threads} cp $fast5_dir/% {output.out_dir}
#         """

# rule sampling_ptr_108_bio3:
#     threads: 10
#     input:
#         filelist = lambda wildcards: query(id="ptr_108_bio3")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/ptr_108_bio3.txt"),
#         out_dir = directory("outputs/sampled_fast5/ptr_108_bio3/")
#     log:
#         "logs/0_read_sampling/sampling_ptr_108_bio3.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         shuf -n 100000 {input.filelist} > {output.sampled_names}
#         mkdir -p {output.out_dir}
#         fast5_dir=$(dirname {input.filelist})
#         cat {output.sampled_names} \
#             | xargs -I % -P {threads} cp $fast5_dir/% {output.out_dir}
#         """

# rule sampling_egr_108_bio1:
#     threads: 10
#     input:
#         filelist = lambda wildcards: query(id="egr_108_bio1")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/egr_108_bio1.txt"),
#         out_dir = directory("outputs/sampled_fast5/egr_108_bio1/")
#     log:
#         "logs/0_read_sampling/sampling_egr_108_bio1.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         shuf -n 100000 {input.filelist} > {output.sampled_names}
#         mkdir -p {output.out_dir}
#         fast5_dir=$(dirname {input.filelist})
#         cat {output.sampled_names} \
#             | xargs -I % -P {threads} cp $fast5_dir/% {output.out_dir}
#         """

# rule sampling_egr_108_bio2:
#     threads: 10
#     input:
#         filelist = lambda wildcards: query(id="egr_108_bio2")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/egr_108_bio2.txt"),
#         out_dir = directory("outputs/sampled_fast5/egr_108_bio2/")
#     log:
#         "logs/0_read_sampling/sampling_egr_108_bio2.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         shuf -n 100000 {input.filelist} > {output.sampled_names}
#         mkdir -p {output.out_dir}
#         fast5_dir=$(dirname {input.filelist})
#         cat {output.sampled_names} \
#             | xargs -I % -P {threads} cp $fast5_dir/% {output.out_dir}
#         """

# rule sampling_egr_108_bio3:
#     threads: 10
#     input:
#         filelist = lambda wildcards: query(id="egr_108_bio3")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/egr_108_bio3.txt"),
#         out_dir = directory("outputs/sampled_fast5/egr_108_bio3/")
#     log:
#         "logs/0_read_sampling/sampling_egr_108_bio3.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         shuf -n 100000 {input.filelist} > {output.sampled_names}
#         mkdir -p {output.out_dir}
#         fast5_dir=$(dirname {input.filelist})
#         cat {output.sampled_names} \
#             | xargs -I % -P {threads} cp $fast5_dir/% {output.out_dir}
#         """

# rule sampling_ptr_109_bio1:
#     threads: 8
#     input:
#         fast5_dir = lambda wildcards: query(id="ptr_109_bio1")["fast5_dir"],
#         filelist = lambda wildcards: query(id="ptr_109_bio1")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/ptr_109_bio1.txt"),
#         out_dir = directory("outputs/sampled_fast5/ptr_109_bio1/")
#     log:
#         "logs/0_read_sampling/sampling_ptr_109_bio1.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         head -n 1 {input.filelist} > {output.sampled_names}
#         shuf -n 100000 {input.filelist} >> {output.sampled_names}
#         mkdir -p {output.out_dir}
#         docker run \
#             {docker_mount_opt} \
#             ccc/ont-fast5-api:latest \
#             fast5_subset \
#                 -i {input.fast5_dir} \
#                 -s {output.out_dir} \
#                 -l {output.sampled_names} \
#                 -n 100000 \
#                 -r \
#                 -t {threads} \
#             2>> {log} \
#             1>> {log}
#         """

# rule sampling_ptr_109_bio2:
#     threads: 8
#     input:
#         fast5_dir = lambda wildcards: query(id="ptr_109_bio2")["fast5_dir"],
#         filelist = lambda wildcards: query(id="ptr_109_bio2")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/ptr_109_bio2.txt"),
#         out_dir = directory("outputs/sampled_fast5/ptr_109_bio2/")
#     log:
#         "logs/0_read_sampling/sampling_ptr_109_bio2.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         head -n 1 {input.filelist} > {output.sampled_names}
#         shuf -n 100000 {input.filelist} >> {output.sampled_names}
#         mkdir -p {output.out_dir}
#         docker run \
#             {docker_mount_opt} \
#             ccc/ont-fast5-api:latest \
#             fast5_subset \
#                 -i {input.fast5_dir} \
#                 -s {output.out_dir} \
#                 -l {output.sampled_names} \
#                 -n 100000 \
#                 -r \
#                 -t {threads} \
#             2>> {log} \
#             1>> {log}
#         """

# rule sampling_egr_109_bio1:
#     threads: 8
#     input:
#         fast5_dir = lambda wildcards: query(id="egr_109_bio1")["fast5_dir"],
#         filelist = lambda wildcards: query(id="egr_109_bio1")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/egr_109_bio1.txt"),
#         out_dir = directory("outputs/sampled_fast5/egr_109_bio1/")
#     log:
#         "logs/0_read_sampling/sampling_egr_109_bio1.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         head -n 1 {input.filelist} > {output.sampled_names}
#         shuf -n 100000 {input.filelist} >> {output.sampled_names}
#         mkdir -p {output.out_dir}
#         docker run \
#             {docker_mount_opt} \
#             ccc/ont-fast5-api:latest \
#             fast5_subset \
#                 -i {input.fast5_dir} \
#                 -s {output.out_dir} \
#                 -l {output.sampled_names} \
#                 -n 100000 \
#                 -r \
#                 -t {threads} \
#             2>> {log} \
#             1>> {log}
#         """

# rule sampling_egr_109_bio2:
#     threads: 8
#     input:
#         fast5_dir = lambda wildcards: query(id="egr_109_bio2")["fast5_dir"],
#         filelist = lambda wildcards: query(id="egr_109_bio2")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/egr_109_bio2.txt"),
#         out_dir = directory("outputs/sampled_fast5/egr_109_bio2/")
#     log:
#         "logs/0_read_sampling/sampling_egr_109_bio2.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         head -n 1 {input.filelist} > {output.sampled_names}
#         shuf -n 100000 {input.filelist} >> {output.sampled_names}
#         mkdir -p {output.out_dir}
#         docker run \
#             {docker_mount_opt} \
#             ccc/ont-fast5-api:latest \
#             fast5_subset \
#                 -i {input.fast5_dir} \
#                 -s {output.out_dir} \
#                 -l {output.sampled_names} \
#                 -n 100000 \
#                 -r \
#                 -t {threads} \
#             2>> {log} \
#             1>> {log}
#         """

# rule sampling_lch_109_bio1:
#     threads: 8
#     input:
#         fast5_dir = lambda wildcards: query(id="lch_109_bio1")["fast5_dir"],
#         filelist = lambda wildcards: query(id="lch_109_bio1")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/lch_109_bio1.txt"),
#         out_dir = directory("outputs/sampled_fast5/lch_109_bio1/")
#     log:
#         "logs/0_read_sampling/sampling_lch_109_bio1.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         head -n 1 {input.filelist} > {output.sampled_names}
#         shuf -n 100000 {input.filelist} >> {output.sampled_names}
#         mkdir -p {output.out_dir}
#         docker run \
#             {docker_mount_opt} \
#             ccc/ont-fast5-api:latest \
#             fast5_subset \
#                 -i {input.fast5_dir} \
#                 -s {output.out_dir} \
#                 -l {output.sampled_names} \
#                 -n 100000 \
#                 -r \
#                 -t {threads} \
#             2>> {log} \
#             1>> {log}
#         """

# rule sampling_lch_109_bio2:
#     threads: 8
#     input:
#         fast5_dir = lambda wildcards: query(id="lch_109_bio2")["fast5_dir"],
#         filelist = lambda wildcards: query(id="lch_109_bio2")["filelist"]
#     output:
#         sampled_names = temp("outputs/sampled_fast5/lch_109_bio2.txt"),
#         out_dir = directory("outputs/sampled_fast5/lch_109_bio2/")
#     log:
#         "logs/0_read_sampling/sampling_lch_109_bio2.log"
#     shell:
#         """
#         echo $(date +%Y%m%d-%H:%M:%S) > {log}
#         head -n 1 {input.filelist} > {output.sampled_names}
#         shuf -n 100000 {input.filelist} >> {output.sampled_names}
#         mkdir -p {output.out_dir}
#         docker run \
#             {docker_mount_opt} \
#             ccc/ont-fast5-api:latest \
#             fast5_subset \
#                 -i {input.fast5_dir} \
#                 -s {output.out_dir} \
#                 -l {output.sampled_names} \
#                 -n 100000 \
#                 -r \
#                 -t {threads} \
#             2>> {log} \
#             1>> {log}
#         """

# rule sampling_ptr_111_bio1:
#     input:
#     output:
#     log:
#     shell:
#         """
#         """

# rule sampling_ptr_111_bio2:
#     input:
#     output:
#     log:
#     shell:
#         """
#         """
