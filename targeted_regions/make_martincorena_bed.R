library(GenomicRanges)

# insert your own GTF path
gencode_file = "gencode.v37lift37.basic.annotation.gtf.gz"
gencode_genes = rtracklayer::import(gencode_file)

# List from Martincorena et al., with the following gene name changes (due to changes in HGNC identifiers over time)
# BAI3 -> ADGRB3, MLL -> KMT2A, MLL2 -> KMT2D, MLL3 -> KMT2C
martincorena_genes = c("ADAM29", "ADAMTS18", "AJUBA", "AKT1", "AKT2", "APOB", "ARID1A", "ARID2", "AURKA", 
              "ADGRB3", "BRAF", "CASP8", "CCND1", "CDH1", "CDKN2A", "CR2", "CREBBP", "CUL3", "DICER1", "EGFR",
              "EPHA2", "ERBB2", "ERBB3", "ERBB4", "EZH2", "FAT1", "FAT4", "FBXW7", "FGFR1", "FGFR2", "FGFR3", "FLG2", 
              "GRIN2A", "GRM3", "HRAS", "IRF6", "KCNH5", "KEAP1", "KRAS", "MET", "KMT2A", "KMT2D", "KMT2C", "MUC17", "NF1", 
              "NFE2L2", "NOTCH1", "NOTCH2", "NOTCH3", "NOTCH4", "NRAS", "NSD1", "PCED1B", "PIK3CA", "PLCB1", "PPP1R3A", 
              "PREX2", "PTCH1", "PTEN", "PTPRT", "RB1", "RBM10", "SALL1", "SCN11A", "SCN1A", "SETD2", "SMAD4", "SMO", 
              "SOX2", "SPHKAP", "SUFU", "TP53", "TP63", "TRIOBP")

# Verify all genes present (they are now)
all(martincorena_genes %in% gencode_genes$gene_name)

# if needed, change from UTR/exon to just exon or something else (you have to decide based on the data source)
our_intervals = gencode_genes[gencode_genes$type %in% c("UTR", "exon") & gencode_genes$gene_name %in% our_genes]

gr = granges(our_intervals)
gr = reduce(sort(unstrand(gr)))
seqlevelsStyle(gr) = 'NCBI'

rtracklayer::export.bed(gr, con = "targeted_regions/martincorena_covered_regions.bed")
