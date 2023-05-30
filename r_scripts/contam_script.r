#! /usr/bin/env Rscript 

library(tidyverse)
library(qiime2R)
library(phyloseq)
library(decontam)
library(biomformat)
library(readr)

#load data into memory
#table_dada2 <- "results/dada2/ASV_table.tsv"
# change table to be consistant with downstream processes
table_dada2 <- "feature-table.tsv"
rooted_tree <- "results/qiime2/phylogenetic_tree/rooted-tree.qza"
taxonomy_file <- "results/qiime2/input/taxonomy.qza"
metadata_file <- "metadata.tsv"

control_name <- read.table("nc_name.txt")

#format metadata in phyloseq compatable
metadata<- read_q2metadata(metadata_file)
metadata <- metadata %>% remove_rownames %>% column_to_rownames(var="SampleID")
metadata_ps <- sample_data(metadata)

# format ASV table into phyloseq compatable
table <- data.frame(read.table(table_dada2,sep='\t'))
#Depricated when using formatted table
#rownames(table) <- table$"X.OTU.ID"
#able <- table %>% select(-c("X.OTU.ID"))
table_ps <- otu_table(table, taxa_are_rows=T)

# format the taxonomy data into phylsoeq compatable
tax_obj <- read_qza(taxonomy_file)
tax_ps <- tax_table(as.matrix(parse_taxonomy(tax_obj$data)))

# format the tree into phyloseq
tree_obj <- read_qza(rooted_tree)
tree <- tree_obj$data

# generate the phyloseq object
ps <- phyloseq(table_ps, tax_ps, tree, metadata_ps)

# add a variable for negative control samples, update list if needed
#sample_data(ps)$is.neg <- sample_data(ps)$Treatment %in% c("NC","Control","Negative Control","control","negative control")
sample_data(ps)$is.neg <- sample_data(ps)$Treatment %in% control_name

# identify the possible contaminats based on 0.5 
# In the prevalence test there is a special value worth knowing, threshold=0.5, 
#that will identify as contaminants all sequences thare are more prevalent in negative controls than in positive samples.

contamdf.prev05 <- isContaminant(ps, method="prevalence", neg="is.neg", threshold=0.5)

#TODO save these contams to disk so we can examine them

ps.contam <- phyloseq::prune_taxa(contamdf.prev05$contaminant,ps)

# filter the phyloseq object to remove contamination
ps.noncontam <- phyloseq::prune_taxa(!contamdf.prev05$contaminant,ps)

# save ASV table to disk

write.table(otu_table(ps.noncontam),sep='\t', "filtered-table.tsv")

# turn phyloseq object into a biom file
biom <- make_biom(data=otu_table(ps.noncontam))

# save biom file to disk
write_biom(biom, "filtered-table.biom")

# save contaminents to disk
sample_data <- psmelt(ps.contam)
otu_table <- otu_table(ps.contam)

contam_together <- merge(otu_table,sample_data)

write.table(contam_together, sep='\t', "contam-features.tsv")
