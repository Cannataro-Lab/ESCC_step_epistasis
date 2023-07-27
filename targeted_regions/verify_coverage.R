library(cancereffectsizeR)
library(data.table)


get_stats = function(dt) {
  return(dt[, .(in_coverage = mean(dist_to_coverage_intervals == 0),
                within_100 = mean(dist_to_coverage_intervals < 101),
                within_1000 = mean(dist_to_coverage_intervals < 1001))])
}

icgc_file = "source_data/raw_data/escc_ucla_2014/data_mutations.txt"
icgc = preload_maf(icgc_file, refset = 'ces.refset.hg19', detect_hidden_mnv = F, 
                   coverage_intervals_to_check = 'targeted_regions/NimbleGen_EZ2_hg19.bed')
icgc_stats = get_stats(icgc)
# Nimblegen EZ 44 Mb panel looks correct; use 100bp padding.


# TCGA ESCA: Nimblegen SeqCap EZ Human Exome Library v3.0
# https://www.nature.com/articles/nature20805
tcga_file = "source_data/tcga_escc.maf" # file is filtered to just ESCC samples
# using a local chain file
tcga = preload_maf(tcga_file, chain_file = '~/reference/chains/hg38ToHg19.over.chain', 
                   refset = 'ces.refset.hg19', detect_hidden_mnv = F, 
                   coverage_intervals_to_check = 'targeted_regions/SeqCap_EZ_Exome_v3_target_regions.bed')

tcga_stats = get_stats(tcga)
# Data has been trimmed to covered regions (not by us; confirmed that original ESCA project MAF is trimmed, too).
# So, padding = 0.

ucla_file = 'source_data/raw_data/escc_ucla_2014/data_mutations.txt'
ucla = preload_maf(ucla_file, refset = 'ces.refset.hg19', detect_hidden_mnv = F,
                   coverage_intervals_to_check = 'targeted_regions/covered_regions_from_UCLA_Lin_S14.bed')
ucla_stats = get_stats(ucla)
# Set padding to 100. (96.3% within 100, vs. 95.7% within coverage.)


# Yuan: Agilent SureSelect Human All Exon V5 (https://academic.oup.com/carcin/article/40/12/1445/5579375)
yuan_file = fread('source_data/yuan.maf')
yuan_maf = preload_maf(maf = yuan_file, refset = 'ces.refset.hg19',
                       detect_hidden_mnv = F, coverage_intervals_to_check = 'targeted_regions/SureSelect_All_Exon_V5_S04380110_Covered_hg19.bed')
yuan_maf = yuan_maf[is.na(problem)] # ignore duplicate records (should be fixed but not a problem for this)
yuan_stats = yuan_maf[, .(in_coverage = mean(dist_to_coverage_intervals == 0),
         within_100 = mean(dist_to_coverage_intervals < 101),
         within_1000 = mean(dist_to_coverage_intervals < 1001))]
# Looks good; use 100bp padding.

# Martincorena (using original data file rather than prepped version)
martincorena_file = as.data.table(readxl::read_excel('source_data/raw_data/martincorena_2018_eso_w_age/aau3879_tables2.xlsx', skip = 16))
martincorena_maf = martincorena_file[, .(Tumor_Sample_Barcode = sampleID, Chromosome = chr, Start_Position = pos, 
                      Reference_Allele = ref, Tumor_Allele = mut)]
martincorena = preload_maf(maf = martincorena_maf, refset = 'ces.refset.hg19',
                           detect_hidden_mnv = F, coverage_intervals_to_check = 'targeted_regions/martincorena_covered_regions.bed')
martincorena = martincorena[is.na(problem)] # some duplicate records; also, some deletions have wrong reference alleles

martincorena_stats = get_stats(martincorena)

# Set padding to 100bp. Note that same number of calls are within 100 and 1000bp. Probably
# Martincorena trimmed to 100bp, and remaining calls are from undisclosed additional sequenced sites.
# In supplementary methods, they say that they also targeted "1,124 SNPs within or around the 74 target genes for 
# targeted copy number analysis." We won't attempt to figure out what they covered outside the 74 genes.


# Yokoyama
yokoyama_file = fread('source_data/yokoyama.maf')
yokoyama = preload_maf(maf = yokoyama_file, refset = 'ces.refset.hg19',
                       detect_hidden_mnv = F, coverage_intervals_to_check = 'targeted_regions/yokoyama_covered_regions.bed')
yokoyama_stats = get_stats(yokoyama)
# 83.0 covered, 88.5% within 100 (set padding 100bp)
# Similar situation to Martincorena: There may be some additional undisclosed targets,
# but we won't do more to try to figure it out.

all_stats = rbindlist(list(tcga = tcga_stats, icgc = icgc_stats, ucla = ucla_stats,
               yuan = yuan_stats, martincorena = martincorena_stats, 
               yokoyama = yokoyama_stats), idcol = 'source')

fwrite(all_stats, 'targeted_regions/variant_call_coverage_stats.txt', sep = "\t")


