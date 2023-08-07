# ESCC_stage_epistasis
## Steps to run the analysis

1.  Obtain data from links described in [input_data/README.md](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/input_data/README.md)

2.  Prep data frames with

    -   [input_data/prep_Liu_maf.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/input_data/prep_Liu_maf.R)
    -   [input_data/prep_Martincorena_maf.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/input_data/prep_Martincorena_maf.R)
    -   [input_data/prep_UCLA_maf.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/input_data/prep_UCLA_maf.R)
    -   [input_data/prep_Yokoyama_maf.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/input_data/prep_Yokoyama_maf.R)
    -   [input_data/prep_Yuan_maf.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/input_data/prep_Yuan_maf.R)
    -   Also be sure to download the ESCC-META dataset from the synapse data portal. More detailed instructions are available within `input_data/README.md`

3.  Recreate cancer effect size analysis by running [analysis/run_analysis.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/analysis/run_analysis.R)

4.  Figures:

    -   Figure 1A: [output_data/Figure_1A_Mutrate.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/output_data/Figure_1A_Mutrate.R)
    -   Figure 1B: [output_data/Figure_1B_cosmic_signatures.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/output_data/Figure_1B_cosmic_signatures.R)
    -   Figure 2: [output_data/Figure_2_selection_coefficient_plot.R](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/output_data/Figure_2_selection_coefficient_plot.R)
    -   Figure 3: [output_data/Figure_3_epistasis_plots](https://github.com/Cannataro-Lab/ESCC_stage_epistasis/blob/main/output_data/Figure_3_epistasis_plots.R)
