# Intro to Snakemake

## 🛠 Prerequisites

Before you begin, make sure you have the following tools installed:

- **Visual Studio Code (VSCode)**  
  Code editor for development and running scripts.  
  👉 [Download VSCode](https://code.visualstudio.com/)

- **Conda or Mamba**  
  Package and environment manager to install required bioinformatics tools.  
  👉 [Download Miniconda](https://docs.conda.io/en/latest/miniconda.html)  
  👉 [Download Mamba (conda faster alternative)](https://mamba.readthedocs.io/en/latest/)

- **WSL (Windows Subsystem for Linux) [Windows only]**  
  Enables Linux tools and environments on Windows.  
  👉 [Set up WSL](https://learn.microsoft.com/en-us/windows/wsl/install)

---

## 📂 Input Files

You will need the following files before running the pipeline:

- **Three FASTQ files**  
  H3K4me3 replicate 1: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/001/SRR1275461/SRR1275461.fastq.gz
  H3K4me3 replicate 2: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/002/SRR1275462/SRR1275462.fastq.gz
  control DNA: ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR127/006/SRR1275466/SRR1275466.fastq.gz

- **Genome reference**  
  [hg19](https://genome-idx.s3.amazonaws.com/bt/hg19.zip)

- **Tools**
  Snakemake
  wget
  fastp
  bowtie2
  macs3
  bedtools2
  multiqc

Make sure these files are available in your working directory or properly linked in your script configuration.

---

