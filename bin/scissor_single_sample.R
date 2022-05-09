#!/usr/bin/env Rscript
"
scissor_single_sample.R

Usage:
    scissor_single_sample.R --bulk_tpm=<bulk_tpm> --sce=<sce> --metadata=<metadata> [options]

Mandatory options:
    --bulk_tpm=<bulk_tpm>       Bulk gene expression as matrix with TPM values (not log transformed!) in rds format
    --sce=<sce>                 SingleCellExperiment in rds format. Must contain raw UMI counts in `X`.
    --metadata=<metadata>       Samplesheet with metadata in TSV format. Binary variables need to be encoded as 0 and 1. Non 0/1 values will be discarded.

The program will take the intersection of the <sample_col> in <metadata> and the colnames of <bulk_tpm>.

Optional options:
    --column=<column>               Column with binary variable for binomial regression
    --tumor_type=<tumor_type>       Restrict analysis to certain tumor type [default: any]
    --surv_time=<surv_time>         Column with survival time for coxph regression. Must be combined with --surv_status
    --surv_status=<surv_status>     Column with survival status for coxpy regression.
    --sample_col=<sample_col>       Column in <metadata> with sample identifier (must match colnames of <bulk_tpm>) [default: sample_id]
    --prefix=<prefix>               Prefix for the output filename
" -> doc

library(conflicted)
library(docopt)

arguments <- docopt(doc, version = "0.1")
print(arguments)

library(Scissor)
library(readxl)
library(SingleCellExperiment)
library(readxl)
library(dplyr)
conflict_prefer("filter", "dplyr")
library(readr)
library(stringr)
library(ggplot2)

sce <- readRDS(arguments$sce)
bulk_tpm <- readRDS(arguments$bulk_tpm)
metadata <- read_tsv(arguments$metadata)
sample_col <- arguments$sample_col
column <- arguments$column
surv_time <- arguments$surv_time
surv_status <- arguments$surv_status
tumor_type <- arguments$tumor_type
prefix <- arguments$prefix

# # For testing only
# sce = readRDS("../data/30_downstream_analyses/scissor/adata_by_patient/full_atlas_merged_lambrechts_thienpont_2018_6149v1_2.rds")
# bulk_tpm = readRDS("../data/14_ici_treatment/Genentech_for_scissor/genentech.rds")
# metadata = read_tsv("../data/14_ici_treatment/Genentech_for_scissor/genentech_clinical_data.tsv")
# sample_col = "sample_id"
# column = NULL
# surv_time = "time_chemo"
# surv_status = "status_chemo"
# column = "response_binary"
# tumor_type = "any"
# surv_time = NULL
# surv_status = NULL

message("subsetting data to intersection")
if (tumor_type != "any") {
    metadata <- metadata %>% filter(type == !!tumor_type)
}
message(paste(dim(metadata), collapse = " "))
common_patients <- sort(intersect(colnames(bulk_tpm), metadata[[sample_col]]))
metadata <- metadata %>%
    filter(!!as.name(sample_col) %in% common_patients) %>%
    arrange(!!as.name(sample_col))
bulk_tpm <- bulk_tpm[, common_patients]
stopifnot(all(colnames(bulk_tpm) == metadata[[sample_col]]))
message(paste(dim(bulk_tpm), collapse = " "))

message("preprocessing single-cell data")
# not ideal to preprocess the same dataset multiple times when testing multiple variables.
# But compared to the scissor regression, the runtime of this part is negligible.
sc_dataset <- Seurat_preprocessing(assays(sce)$X)
p <- DimPlot(sc_dataset, reduction = "umap", label = T, label.size = 5)
ggsave("dim_plot_seurat.pdf", plot = p)

if (!is.null(column)) {
    message("running binomial regression")
    message("subsetting to samples that have values of 0 or 1 in the selected clinical variable")
    sample_prefix <- column
    metadata_subset <- metadata %>% filter(!!as.name(column) %in% c(0, 1))
    bulk_tpm_subset <- bulk_tpm[, metadata_subset[[sample_col]]]
    stopifnot(all(colnames(bulk_tpm_subset) == metadata_subset[[sample_col]]))
    message(paste(dim(bulk_tpm_subset), collapse = " "))

    infos1 <- Scissor(bulk_tpm_subset, sc_dataset, metadata_subset[[column]],
        tag = c(0, 1),
        alpha = sqrt(2)^-(24:2),
        cutoff = 0.3,
        family = "binomial", Save_file = "scissor.RData"
    )

    Scissor_select <- rep(NA, ncol(sc_dataset))
    names(Scissor_select) <- colnames(sc_dataset)

    # sanity check due to some weirdness we observed
    stopifnot("too many scissor cells" = length(infos1$Scissor_neg) <= dim(sce)[2])
    stopifnot("too many scissor cells" = length(infos1$Scissor_pos) <= dim(sce)[2])

    Scissor_select[infos1$Scissor_pos] <- "scissor+"
    Scissor_select[infos1$Scissor_neg] <- "scissor-"
} else if (!is.null(surv_time) && !is.null(surv_status)) {
    message("running coxph regression")
    metadata_subset <- metadata %>% filter(!is.na(!!as.name(surv_time)), !is.na(!!as.name(surv_status)))
    bulk_tpm_subset <- bulk_tpm[, metadata_subset[[sample_col]]]
    stopifnot(all(colnames(bulk_tpm_subset) == metadata_subset[[sample_col]]))
    message(paste(dim(bulk_tpm_subset), collapse = " "))

    phenotype <- metadata_subset %>% select(time = !!as.name(surv_time), status = !!as.name(surv_status))
    sample_prefix <- paste0(surv_status, "_", surv_time)
    infos1 <- Scissor(bulk_tpm_subset, sc_dataset, phenotype,
        alpha = sqrt(2)^-(24:2),
        cutoff = 0.3,
        family = "cox", Save_file = "scissor.RData"
    )

    Scissor_select <- rep(NA, ncol(sc_dataset))
    names(Scissor_select) <- colnames(sc_dataset)

    # sanity check due to some weirdness we observed
    stopifnot("too many scissor cells" = length(infos1$Scissor_neg) <= dim(sce)[2])
    stopifnot("too many scissor cells" = length(infos1$Scissor_pos) <= dim(sce)[2])

    # it is a bit counterintuitive, but scissor+ -> worse survival, scissor- -> better survival
    Scissor_select[infos1$Scissor_pos] <- "worse survival"
    Scissor_select[infos1$Scissor_neg] <- "better survival"
} else {
    stop("Either specify --column or both --surv_time and --surv_status!")
}


message("writing output")
sc_dataset_meta <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
p <- DimPlot(sc_dataset_meta,
    reduction = "umap", group.by = "scissor",
    cols = c("royalblue", "indianred1"), pt.size = .2, order = c(2, 1)
)
ggsave("dim_plot_scissor.pdf", plot = p)

write_tsv(as.data.frame(Scissor_select) %>% as_tibble(rownames = "cell_id"), sprintf("scissor_%s_%s_%s.tsv", prefix, tumor_type, sample_prefix))
