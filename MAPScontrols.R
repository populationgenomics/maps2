# LT 21/10/21

# Rewrote the hail extraction of MAPS calibration data
# hail script MAPSextract.py 
# generated file syn_ps.tsv

# applied it to gnomAD controls (n=)
# hail script MAPSextract_controls.py 
# generated file syn_ps_controls.tsv

# 28/10/2021
# add non cancer and other subsets


# compare with MAPS fit 10/2019
require(tidyverse)

syn_orig = read.table("../syn_ps.filtered.singleton.tsv",header=TRUE)
syn = read.table("syn_ps.tsv",header=TRUE)

plot(ps~mu, data =syn, xlab="mutability", ylab="proportion of singletons")
points(ps~mu, data =syn_orig, col=2)
# same

# confirmed
plot(syn_orig$ps ~syn$ps, main='Check new MAPS code', xlab='Prop. singletons (original)', ylab='Prop. singletons (new code)')
abline(a=0,b=1)

# now load the controls
syn_ctrl = read.table("syn_ps_controls.tsv",header=TRUE)

plot(ps~mu, data =syn, log='',xlab="mutability", ylab="proportion of singletons")
points(ps~mu, data =syn_ctrl, col=2)
legend('topright', legend=c("gnomAD exomes (n=141,456)","controls (n=60,146)"), col=c(1,2), lwd=2, bty='n')

# add genomes
syn_gen = read.table("../syn_ps.genomes.tsv",header=TRUE)
#points(ps~mu, data =syn_gen, col=3)
#legend('topright', legend=c("gnomAD exomes (n=141,456)","controls (n=60,146)", "v2.1.1 genomes (n=)"), col=c(1,2,3), lwd=2, bty='n')