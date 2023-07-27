library(GenomicRanges)

# insert your own GTF path
gencode_file = "gencode.v37lift37.basic.annotation.gtf.gz"
gencode_genes = rtracklayer::import(gencode_file)

# Listed in table S7 (see source_data/raw_data/yokoyama/41586_2018_811_MOESM3_ESM.xlsx)
yokoyama_genes = c("NOTCH1", "TP53", "PIK3CA", "PPM1D", "FAT1", "ZFP36L2", "EP300", "NOTCH2", "NOTCH3",
                  "CHEK2", "KMT2D", "PAX9", "ZNF750", "CUL3", "CREBBP", "AJUBA", "NFE2L2", "PLXNB2",
                  "PTCH1", "CDKN2A", "FBXW7", "KDM6A", "RB1", "TGFBR2")

all(yokoyama_genes %in% gencode_genes$gene_name)

# if needed, change from UTR/exon to just exon or something else (you have to decide based on the data source)
our_intervals = gencode_genes[gencode_genes$type %in% c("UTR", "exon") & gencode_genes$gene_name %in% our_genes]

gr = granges(our_intervals)
gr = reduce(sort(unstrand(gr)))
seqlevelsStyle(gr) = 'NCBI'

rtracklayer::export.bed(gr, con = "targeted_regions/yokoyama_covered_regions.bed")

