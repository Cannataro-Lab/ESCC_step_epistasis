================================================================================
Single-Cell RNA-seq Analysis Pipeline for GSE160269
================================================================================

STEP 1: Get the data from GEO
------------------------------
You'll need to download these files from GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE160269):

Required files:
  - GSE160269_series_matrix.txt (metadata)
  - GSE160269_CD45neg_UMIs.txt.gz (UMI counts for CD45 negative cells)
  - GSE160269_CD45pos_UMIs.txt.gz (UMI counts for CD45 positive cells)

Optional files (referenced but not strictly required):
  - GSE160269_CD45neg_cells.txt.gz
  - GSE160269_CD45pos_cells.txt.gz

Put all these files in the same directory as the analysis scripts.


STEP 2: Set up your environment
--------------------------------
Build and load the conda environment using the environment.yml file:

  conda env create -f environment.yml
  conda activate escc_analysis

This will install all the necessary packages (scanpy, pandas, numpy, etc.)


STEP 3: Combine the UMI files
------------------------------
The data comes split into CD45+ and CD45- populations but we aren't interested in that question, so first you need to merge them.

Run the combining script:

  python mem_eff_combine.py

This will create a file called 'all_combined_UMIs.txt.gz' with everything merged together which will be the basis of the rest of the analysis. 

If you want to change the input/output filenames, check out the main function at the bottom of mem_eff_combine.py where you can modify the file paths.


STEP 4: Run the main analysis
------------------------------
Now we can run the actual single-cell analysis:

  python scrna_analysis_parallel.py

This script does QC, normalization, clustering, UMAP, differential expression, and a basic pseudotime analysis. It'll save a post-processed matrix after all the cleaning steps, which means if you want to tweak the downstream analyses later, subsequent runs will be much, much faster since it'll skip the preprocessing.

A few things you can customize (at the top of the script):
  - N_THREADS
  - Verbosity and other parameters

But some important things must get modifyined in the main function itself:
  -  input/output file paths
  - Modify the genes of interest to highlight in plots (just edit the python list)


STEP 5: (Optional) More sophisticated pseudotime analysis
----------------------------------------------------------
The main script does a basic pseudotime analysis, but if you want it done by cell type in a more nuanced way, run:

  python pseudotime_reanalysis.py

This assumes you've already run scrna_analysis_parallel.py and have the post-processed data file in your directory. It'll do trajectory analysis with different root anchoring strategies for each cell type.
The scripts expect files in the same directory by default, so keep everything together unless you modify the paths 