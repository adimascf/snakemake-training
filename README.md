# ChIP-seq Snakemake Pipeline

This repository contains a Snakemake workflow for processing ChIP-seq data, specifically for the H3K4me3 histone mark.  
The pipeline takes raw FASTQ files, processes them through multiple analysis steps, and outputs final bigwig and peak files ready for visualization and downstream analysis.

---

## 🛠 What Does This Pipeline Do?

The Snakemake workflow runs the following steps:

1. **Download raw sequencing data (download_fastq)**  
   Retrieves replicate and control FASTQ files directly from the SRA FTP server.

2. **Quality control and filtering (run_fastp)**  
   Uses `fastp` to trim, filter, and generate quality control reports for each replicate and control.

3. **Download reference genome (download_refs)**  
   Downloads the hg19 Bowtie2 index files required for mapping.

4. **Read mapping (read_mapping)**  
   Aligns each replicate and control FASTQ file to the hg19 genome using `bowtie2`.  
   Filters and sorts the aligned reads using `samtools`.

5. **Peak calling (call_peaks)**  
   Uses `MACS3` to:
   - Call peaks using both replicates together vs. control.
   - Call peaks for each replicate individually vs. control.
   - Generate pileup bedGraph files.

6. **reproducibility_analysis: Reproducibility analysis (reproducibility_analysis)**  
   Uses `idr` (Irreproducible Discovery Rate) to:
   - Assess consistency between replicate peak sets.
   - Filter reliable peaks based on IDR thresholds.

7. **create_bigwig: BigWig generation (create_bigwig)**  
   Converts pileup bedGraph files into bigWig files for genome browser visualization.

---

## 🔧 How to use and run the Workflow

1. **Download and unzip this repository**

2. **Open the project folder in VS Code**

3. **activate snakemake env**
   ```bash
   conda activate snakemake

4. **Run the pipeline:**
   ```bash
   snakemake --cores --use-conda
