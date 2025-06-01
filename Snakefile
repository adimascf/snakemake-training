rule all:
    input:
        "data/processed/peak_calling/h3k4me3.peak.bw"

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
        directory("data/genome/hg19")
    shell:
        """
        mkdir -p data/genome/hg19
        wget -O data/genome/hg19.zip {params.url}
        unzip data/genome/hg19.zip -d data/genome/hg19
        """

rule read_mapping:
    input:
        rep1 = "data/processed/fastp/h3k4me3_rep1.fastq.gz",
        rep2 = "data/processed/fastp/h3k4me3_rep2.fastq.gz",
        control = "data/processed/fastp/control.fastq.gz",
        reference = "data/genome/hg19"
    log:
        rep1 = "data/processed/fastp/h3k4me3_rep1.log",
        rep2 = "data/processed/fastp/h3k4me3_rep2.log",
        control = "data/processed/fastp/control.log"
    conda:
        "envs.yaml"
    output:
        rep1 = "data/processed/mapped/h3k4me3_rep1.bam",
        rep2 = "data/processed/mapped/h3k4me3_rep2.bam",
        control = "data/processed/mapped/control.bam"
    shell:
        """
        mkdir -p data/processed/mapped/
        # untuk replicate 1
        bowtie2 -p 4 -x {input.reference}/hg19 -U {input.rep1} 2> {log.rep1} \
            | samtools view -@ 2 -F 4 -q 1 -bS \
            | samtools sort -@ 2 -o {output.rep1}

         # untuk replicate 2
        bowtie2 -p 4 -x {input.reference}/hg19 -U {input.rep2} 2> {log.rep2} \
            | samtools view -@ 2 -F 4 -q 1 -bS \
            | samtools sort -@ 2 -o {output.rep2}

         # untuk control
        bowtie2 -p 4 -x {input.reference}/hg19 -U {input.control} 2> {log.control} \
            | samtools view -@ 2 -F 4 -q 1 -bS \
            | samtools sort -@ 2 -o {output.control}
        """

rule call_peaks:
    input:
        rep1 = "data/processed/mapped/h3k4me3_rep1.bam",
        rep2 = "data/processed/mapped/h3k4me3_rep2.bam",
        control = "data/processed/mapped/control.bam"
    params:
        geneme_size = 2.7e9,
        outdir = "data/processed/peak_calling/",
        both = "h3k4me3_both",
        rep1 = "h3k4me3_rep1",
        rep2 = "h3k4me3_rep2"
    conda:
        "macs3.yaml"
    output:
        "data/processed/peak_calling/h3k4me3_both_peaks.narrowPeak",
        "data/processed/peak_calling/h3k4me3_rep1_peaks.narrowPeak",
        "data/processed/peak_calling/h3k4me3_rep2_peaks.narrowPeak",
        "data/processed/peak_calling/h3k4me3_rep1_treat_pileup.bdg",
        "data/processed/peak_calling/h3k4me3_rep2_treat_pileup.bdg",
        "data/processed/peak_calling/h3k4me3_both_treat_pileup.bdg"
    shell:
        """
        mkdir -p data/processed/peak_calling/
        # first run, both replicates against control
        macs3 callpeak -t {input.rep1} {input.rep2} -c {input.control} \
            -f BAM -g {params.geneme_size} -n {params.both} -B \
            --outdir {params.outdir}

        # second run, rep 1 aganst control
        macs3 callpeak -t {input.rep1} -c {input.control} \
            -f BAM -g {params.geneme_size} -n {params.rep1} -B \
            --outdir {params.outdir}

        # third run, rep 2 against control
        macs3 callpeak -t {input.rep2} -c {input.control} \
            -f BAM -g {params.geneme_size} -n {params.rep2} -B \
            --outdir {params.outdir}
        """

rule reproducibility_analysis:
    input:
        both = "data/processed/peak_calling/h3k4me3_both_peaks.narrowPeak",
        rep1 = "data/processed/peak_calling/h3k4me3_rep1_peaks.narrowPeak",
        rep2 = "data/processed/peak_calling/h3k4me3_rep2_peaks.narrowPeak"
    params:
        log = "data/processed//peak_calling/h3k4me3.idr.log",
        intermediate_out = "data/processed/peak_calling/h3k4me3_idrValues"
    conda:
        "idr.yaml"
    output:
        "data/processed/peak_calling/h3k4me3.idr.final.bed"
    shell:
        """
        # Run the IDR from macs outputs
        idr --samples {input.rep1} {input.rep2} \
            --peak-list {input.both} \
            --input-file-type narrowPeak \
            -o {params.intermediate_out} \
            --plot --use-best-multisummit-IDR \
            --log-output-file {params.log}

        # Filter out unreliable peaks
        awk -F '\t' '$12 > -log(0.05)/log(10)' {params.intermediate_out} | cut -f1-10 > {output}
        """

rule create_bigwig:
    input:
        both = "data/processed/peak_calling/h3k4me3_both_treat_pileup.bdg"
    params:
        chrom_size = "http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes"
    conda:
        "bigwig.yaml"
    output:
        "data/processed/peak_calling/h3k4me3.peak.bw"
    shell:
        """
        wget -O data/raw/chrom.size {params.chrom_size}
        bedGraphToBigWig {input.both} data/raw/chrom.size {output}
        """
