#! /usr/bin/env Rscript 

library(tidyverse)
library(qiime2R)
library(phyloseq)
library(decontam)
library(biomformat)

#load data into memory
table_dada2 <- "results/dada2/ASV_table.tsv"
rooted_tree <- "results/qiime2/phylogenetic_tree/rooted-tree.qza"
taxonomy_file <- "results/qiime2/input/taxonomy.qza"
metadata_file <- "metadata.tsv"

#create phyloseq Obj
metadata<- read_q2metadata(metadata_file)
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="SampleID")
metadata_ps <- sample_data(metadata)

table <- data.frame(read_table(table_dada2))
rownames(table) <- table$ASV_ID
table <- table %>% select(-c("ASV_ID"))
table_ps <- otu_table(table, taxa_are_rows=T)

tax_obj <- read_qza(taxonomy_file)
tax_ps <- tax_table(as.matrix(parse_taxonomy(tax_obj$data)))


tree_obj <- read_qza(rooted_tree)
tree <- tree_obj$data

ps <- phyloseq(table_ps, tax_ps, tree, metadata_ps)



#TODO filter out Negative control samples

#This is biom v1.0 can we turn it into a hd5 format?
biom <- make_biom(otu_table(ps))

write_biom(biom, "biomtest.biom")
#Save table as a qza
write.table(OTUdf, file='table.tsv', quote=FALSE, sep='\t', row.names=F)