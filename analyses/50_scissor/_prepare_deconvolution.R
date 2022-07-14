#!/usr/bin/env Rscript

library(conflicted)
library(dplyr)
conflict_prefer("filter", "dplyr")
library(tidyr)
library(tibble)
library(readr)
library(immunedeconv)
library(BiocParallel)
conflict_prefer("deconvolute", "immunedeconv")
BiocParallel::register(MulticoreParam(10))

set_cibersort_binary("/home/sturm/projects/2022/CIBERSORT/CIBERSORT.R")
set_cibersort_mat("/home/sturm/projects/2022/CIBERSORT/LM22.txt")
cell_types_of_interest = c(
  "B cell",
  "B cell plasma",
  "Cancer associated fibroblast",
  "Endothelial cell",
  "Macrophage",
  "Mast cell",
  "Monocyte",
  "Myeloid dendritic cell",
  "Neutrophil",
  "NK cell",
  "Plasmacytoid dendritic cell",
  "T cell CD4+",
  "T cell CD8+"
)

tcga_tpm = read_rds("../../data/13_tcga/for_scissor/nsclc_primary_tumor.rds")
clinical_data = read_tsv("../../tables/tcga/clinical_data_for_scissor.tsv")
patient_barcodes = intersect(colnames(tcga_tpm), clinical_data$TCGA_patient_barcode)

clinical_data = clinical_data |> filter(TCGA_patient_barcode %in% patient_barcodes)
tpm_mat = tcga_tpm[, patient_barcodes]

deconv_result = BiocParallel::bplapply(deconvolution_methods, function(method) {
  deconvolute(tpm_mat, method=method, indications = clinical_data$type, expected_cell_types=cell_types_of_interest) |>
    mutate(method = method, .after=cell_type)
}) |> bind_rows()

deconv_result_remapped = lapply(deconvolution_methods, function(method) {
  deconv_result |> filter(method==!!method) |> select(-method) |> map_result_to_celltypes(cell_types_of_interest) |> as_tibble(rownames="cell_type") |> mutate(method=method)
}) |> bind_rows()

write_tsv(deconv_result, "../../data/13_tcga/immunedeconv/nsclc_deconvolution.tsv")
write_tsv(deconv_result_remapped, "../../data/13_tcga/immunedeconv/nsclc_deconvolution_remapped.tsv")

