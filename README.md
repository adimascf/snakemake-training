## 🔧 Setup Instructions

Until `read_mapping` rule

1. **Download and unzip the repository**  
   - Download the ZIP file from this repository.  
   - Unzip the file to your desired location.

2. **Open the folder in VSCode (optional but recommended)**  
   - Launch VSCode and open the unzipped folder.

3. **Activate the Snakemake environment**  
   Open your terminal and run:
   ```bash
   conda activate snakemake
   ```
4. **Run Snakemake** 
   ```bash
   snakemake --cores 4 --use-cond read_mapping
   ```