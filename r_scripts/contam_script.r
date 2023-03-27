#! /usr/bin/env Rscript 

library(tidyverse)
library(qiime2R)
library(phyloseq)
library(decontam)

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
nice_tax <- tax_table(matrix(parse_taxonomy(tax_obj$data)))
tax_ps <- tax_table(as.matrix(nice_tax))

tree_obj <- read_qza(rooted_tree)
tree <- tree_obj$data

ps <- phyloseq(table_ps, tax_ps, tree, metadata_ps)

#TODO filter out Negative control samples

#pre-processing filtered table
OTU1 = as(otu_table(table_phylo), "matrix")
#if(taxa_are_rows(table_phylo)){OTU1 <- t(OTU1)}
# Coerce to data.frame
OTUdf = as.data.frame(OTU1)
OTUdf$"#OTU" <- rownames(OTUdf)
rownames(OTUdf) <- c()

#Save table as a qza
write.table(OTUdf, file='table.tsv', quote=FALSE, sep='\t', row.names=F)