# Data sources used for the analysis

**Only the data sets from Liu et al., Yuan et al. and ESCC-META need to be downloaded separately. Save all of these files to the folder input_data/raw_data/. The other data sets can be pulled by running the data_prep scripts.**

UCLA data set was obtained from <https://cbioportal-datahub.s3.amazonaws.com/escc_ucla_2014.tar.gz> on June 8 2022

Run `/input_data/prep_UCLA_maf.R` to download and prepare these data

Martincorena et al. dataset was obtained from <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6298579/bin/NIHMS80426-supplement-Supplementary_Table_2.xlsx> on June 29 2022


Run `/input_data/prep_Martincorena_maf.R` to download and prepare these data

Yokoyama et al. dataset was obtained from <https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0811-x/MediaObjects/41586_2018_811_MOESM3_ESM.xlsx> on June 12 2022


Run `/input_data/prep_Yokoyama_maf.R` to download and prepare these data

**Yuan et al. dataset was obtained from <https://academic.oup.com/carcin/article/40/12/1445/5579375?login=true> (Supplementary table 3) on July 12 2022**


 Download Yuan et al. dataset from https://academic.oup.com/carcin/article/40/12/1445/5579375?login=true  ( [Supplementary table 3](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/carcin/40/12/10.1093_carcin_bgz162/1/bgz162_suppl_supplementary_table_3.xlsx?Expires=1694441539&Signature=cShag6hxbf7rr5PKFhUi6HFTMu5au5N~X2U-6E13x6lmcEX0Bo8oFERVU6EF4f3COJifR2K3dzu6bLGu4CHiluM1lFmvSWu5tHqW8cdYF1h5EdQjIPNI1PcCET5mAjHgqlJDZIgkljVFoEUkcL8YvztPE2EdznZ3RqVOzHL9wFynbn5BmPxZL~yjoQb1hzRu3xHU52FS9XiVzYK9jVm22cHvl4JSlKAkVp2Kvw40siZwQ0DAobx-w1oMOckzkPIOc2bR~2J2T0QwC-nWeALx7P94NNguoonTyRY34uHC5GylTpqskVpwLbmfYe6njvfZFhjrHC~uSzA39eMP0zA5UQ__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA) )

Save file to input_data/raw_data/


**ESCC-META dataset which included 18 different datasets was obtained from <https://www.synapse.org/#!Synapse:syn33401857> on October 19 2022**


Save these data to `input_data/ESCC-META/` and extract within that directory


**Liu et al. dataset was obtained from <https://www.gastrojournal.org/cms/10.1053/j.gastro.2017.03.033/attachment/37b9d627-f726-495b-8143-71c5f0e73f40/mmc2.xlsx> on June 29 2023**

Download Liu et al. dataset from https://www.gastrojournal.org/cms/10.1053/j.gastro.2017.03.033/attachment/37b9d627-f726-495b-8143-71c5f0e73f40/mmc2.xlsx
Save file to input_data/raw_data/

Then run `input_data/prep_Liu_maf.R`
