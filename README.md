## 🔧 Setup Instructions

Up until `download_refs` rule

1. **Download and unzip the repository**  
   - Download the ZIP file from this repository.  
   - Unzip the file to your desired location.

2. **Open the folder in VSCode**  
   - Launch VSCode and open the unzipped folder.

3. **Activate the Snakemake environment**  
   Open your terminal and run:
   ```bash
   conda activate snakemake
   ```
4. **Run Snakemake** 
   ```bash
   snakemake --cores 4 --use-cond download_refs
   ```
⏱ Note: Downloading human reference genome (~3GB) may take approximately 1 hour to complete, depending on your system’s performance and available resources.
