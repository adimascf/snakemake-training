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
        hg19 = "data/genome/hg19"
    shell:
        """
        mkdir -p data/genome/hg19
        wget -O data/genome/hg19.zip {params.url}
        unzip data/genome/hg19.zip -d data/genome/hg19
        """
