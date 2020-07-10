#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager")
  }

BiocManager::install("simpleaffy")
BiocManager::install("mouse4302.db")
# load data
print(args)
library(simpleaffy)
library(mouse4302.db)
data <- read.affy("pheno.txt", args[1])
data_rma <- call.exprs(data, "rma")

data_rma_df <- data.frame(exprs(data_rma))
annot <- data.frame(
    ACCNUM = sapply(contents(mouse4302ACCNUM), paste, collapse = ", "),
    SYMBOL = sapply(contents(mouse4302SYMBOL), paste, collapse = ", "),
    DESC = sapply(contents(mouse4302GENENAME), paste, collapse = ", ")
)

all <- merge(annot, data_rma_df, by.x = 0, by.y = 0, all = T)
write.table(all, file = "exprs.ann.txt", sep = "\t")