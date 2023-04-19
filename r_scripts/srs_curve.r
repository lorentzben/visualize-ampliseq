#! /usr/bin/env Rscript
library(SRS)
library(qiime2R)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("At least one argument must be supplied", call.=FALSE)
} else if (args[1] != 0) {
  cmax <- args[1]

  table = read_qza('table.qza')
  just_table <- table$data

  png('SRS_curve.png')
  srs_curve = SRScurve(just_table, cmax)
  dev.off()

  Artifact.export_data(srs_curve,'srs')
}
