#!/usr/bin/env Rscript
library(SCEVAN)
library(SingleCellExperiment)

args <- commandArgs(trailingOnly = TRUE)

data_mat <- readRDS(args[1])
cores <- as.integer(args[2])
name <- args[3]
data_mat <- assay(data_mat, "X")
data_mat <- as.matrix(data_mat)


res <- SCEVAN::pipelineCNA(count_mtx = data_mat, sample = name, par_cores = cores, SUBCLONES = FALSE)
write.csv(res, file = "scevan_result.csv")
