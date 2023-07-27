library(data.table)

# File is Supplementary Table 14 from 10.1038/ng.2935 (Lin et al.)
dt = as.data.table(readxl::read_excel('targeted_regions/41588_2014_BFng2935_MOESM32_ESM.xlsx', skip = 1))
dt[, chr := gsub('^chr', '', Chromsome)] # typo in table
gr = makeGRangesFromDataFrame(dt[, .(chr, start = Start, end = End)])
export(reduce(sort(gr)), 'covered_regions_from_UCLA_Lin_S14.bed')
