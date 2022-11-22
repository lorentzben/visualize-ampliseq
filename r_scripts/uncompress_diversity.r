#! /usr/bin/env Rscript

library(qiime2R)

# read in bray curtis distance

bray <-  qiime2R::read_qza("bray_curtis_distance_matrix.qza")

write.table(as.matrix(bray$data),sep='\t', file='bray_distance.tsv', col.names=NA)

# make jaccard div matrix
jaccard <- qiime2R::read_qza("jaccard2.qza")
write.table(as.matrix(jaccard$data),sep='\t', file='jaccard_distance.tsv', col.names=NA)

# make unweighted div matrix
u_uni <- qiime2R::read_qza("unweighted_unifrac_distance_matrix.qza")
write.table(as.matrix(u_uni$data),sep='\t', file='unweighted_unifrac_distance.tsv', col.names=NA)

# make weighted div matrix
w_uni <- qiime2R::read_qza("weighted_unifrac_distance_matrix.qza")
write.table(as.matrix(w_uni$data),sep='\t', file='weighted_unifrac_distance.tsv', col.names=NA)