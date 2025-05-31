rule download_fastq:
    params:
        rep1 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/001/SRR1275461/SRR1275461.fastq.gz",
        rep2 = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/002/SRR1275462/SRR1275462.fastq.gz",
        control = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/006/SRR1275466/SRR1275466.fastq.gz"
    output:
        rep1 = "data/raw/fastqs/h3k4me3_rep1.fastq.gz",
        rep2 = "data/raw/fastqs/h3k4me3_rep2.fastq.gz",
        control = "data/raw/fastqs/control.fastq.gz"
    shell:
        """
        wget -O {output.rep1} {params.rep1}
        wget -O {output.rep2} {params.rep2}
        wget -O {output.control} {params.control}
        """

rule run_fastp:
    input:
        rep1 = "data/raw/fastqs/h3k4me3_rep1.fastq.gz",
        rep2 = "data/raw/fastqs/h3k4me3_rep2.fastq.gz",
        control = "data/raw/fastqs/control.fastq.gz"
    log:
        rep1 = "data/processed/fastp/h3k4me3_rep1",
        rep2 = "data/processed/fastp/h3k4me3_rep2",
        control = "data/processed/fastp/control"
    output:
        rep1 = "data/processed/fastp/h3k4me3_rep1.fastq.gz",
        rep2 = "data/processed/fastp/h3k4me3_rep2.fastq.gz",
        control = "data/processed/fastp/control.fastq.gz"
    conda:
        "envs.yaml"
    shell:
        """
        mkdir -p data/processed/fastp
        fastp --in1 {input.rep1} --out1 {output.rep1} --html {log.rep1}.html --json {log.rep1}.json
        fastp --in1 {input.rep2} --out1 {output.rep2} --html {log.rep2}.html --json {log.rep2}.json
        fastp --in1 {input.control} --out1 {output.control} --html {log.control}.html --json {log.control}.json
        """
rule download_refs:
    params:
        url = "https://genome-idx.s3.amazonaws.com/bt/hg19.zip"
    output:
        hg19 = "data/genome/hg19/"
    shell:
        """
        mkdir -p data/genome/hg19
        wget -O {output.hg19} {params.url}
        unzip {output.hg19} -d data/genome/hg19
        """

rule read_mapping:
    input:
        rep1 = "data/processed/fastp/h3k4me3_rep1.fastq.gz",
        rep2 = "data/processed/fastp/h3k4me3_rep2.fastq.gz",
        control = "data/processed/fastp/control.fastq.gz",
        genome = "data/genome/hg19/"
    threads: 4
    log:
        rep1 = "data/processed/mapped/h3k4me3_rep1.log",
        rep2 = "data/processed/mapped/h3k4me3_rep2.log",
        control = "data/processed/mapped/control.log"
    conda:
        "envs.yaml"
    output:
        rep1 = "data/processed/mapped/h3k4me3_rep1.bam",
        rep2 = "data/processed/mapped/h3k4me3_rep2.bam",
        control = "data/processed/mapped/control.bam"
    shell:
        """
        mkdir -p data/processed/mapped
        bowtie2 -p {threads} -x {input.genome}/hg19 -U {input.rep1} 2> {log.rep1}.log \
            | samtools view -@ 2 -F 4 -q 1 -bS \
            | samtools sort -@ 2 -o {output.rep1}

        bowtie2 -p {threads} -x {input.genome}/hg19 -U {input.rep2} 2> {log.rep2}.log \
            | samtools view -@ 2 -F 4 -q 1 -bS \
            | samtools sort -@ 2 -o {output.rep2}

        bowtie2 -p {threads} -x {input.genome}/hg19 -U {input.control} 2> {log.control}.log \
            | samtools view -@ 2 -F 4 -q 1 -bS \
            | samtools sort -@ 2 -o {output.control}
        """