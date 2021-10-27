# LT 22/10/2021

#
# first look at SFS with mutability
require(tidyverse)

syn = read.table("syn_sfs.tsv",header=TRUE)
syn = filter(syn, protein_coding=='true')

syn= filter(syn, mu == max(syn$mu))

# filter binned categories
syn = filter(syn, cat<128)
plot(variant_count/sum(variant_count) ~cat, data=syn,
     main = 'ACG -> ATG, high methylation', xlab='Allele count', ylab='Proportion')
