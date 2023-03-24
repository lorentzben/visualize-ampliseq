#! /usr/bin/env Rscript 

library(tidyverse)
library(qiime2R)
library(phyloseq)
library(decontam)

#load data into memory
table_dada2 <- "results/qiime2/input/table.qza"
rooted_tree <- "results/qiime2/phylogenetic_tree/rooted-tree.qza"
taxonomy_file <- "results/qiime2/input/taxonomy.qza"
metadata_file <- "metadata.tsv"

#create phyloseq Obj
table_phylo <- qza_to_phyloseq(table_dada2,rooted_tree,taxonomy_file,metadata_file)

#TODO filter out Negative control samples

#pre-processing filtered table
OTU1 = as(otu_table(table_phylo), "matrix")
#if(taxa_are_rows(table_phylo)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
OTUdf$"#OTU ID" <- rownames(OTUdf)
rownames(OTUdf) <- c()

#Save table as a qza
write.table(OTUdf, file='table.tsv', quote=FALSE, sep='\t', row.names=F)