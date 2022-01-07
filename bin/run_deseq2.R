#!/usr/bin/env Rscript

# Import library
library("BiocParallel")
library("DESeq2")
library("org.Hs.eg.db")
library("dplyr")
library("IHW")
library("tibble")
library("readr")
library("argparser", quietly = TRUE)


# Create a parser
p <- arg_parser("Input for Differential Expression analysis")

# Add command line arguments
p <- add_argument(p, "countData", help = "Count Matrix with non-negative integers", type = "character")
p <- add_argument(p, "colData", help = "Matrix with sample annotations", type = "character")
p <- add_argument(p, "--covariate_formula", help = "Additional covariates that will be added to the design formula", default = "")
p <- add_argument(p, "--resDir", help = "Output result directory", default = "./results")
p <- add_argument(p, "--cond_col", help = "Column in sample annotation that contains the condition", default = "group")
p <- add_argument(p, "--paired_grp", help = "Column containing the name of the paired samples, when dealing with paired data", default = "donor")
p <- add_argument(p, "--prefix", help = "Prefix of result file", default = "test")
p <- add_argument(p, "--sample_col", help = "Column in sample annotation that contains the sample names", default = "sample")
p <- add_argument(p, "--id_type", help = "Type of the identifier in the `gene_id` column compatible with AnnotationDbi", default = "ENSEMBL")
p <- add_argument(p, "--fdr", help = "False discovery rate for GO analysis and volcano plots", default = 0.1)
p <- add_argument(p, "--c1", help = "Contrast level 1 (perturbation)", default = "grpA")
p <- add_argument(p, "--c2", help = "Contrast level 2 (baseline)", default = "grpB")
p <- add_argument(p, "--cpus", help = "Number of cpus", default = 8)

# Parse the command line arguments
argv <- parse_args(p)

# Function for removing the ENSG version
remove_ensg_version <- function(x) gsub("\\.[0-9]*$", "", x)

# test data path
# sampleAnnotationCSV = "../tests/testdata/example_nfcore/rnaseq_samplesheet_nfcore-3.1.csv"
# readCountFile = "../tests/testdata/example_nfcore/salmon.merged.gene_counts.subset.tsv"

results_dir <- argv$resDir
paired_grp <- argv$paired_grp
prefix <- argv$prefix
cond_col <- argv$cond_col
sample_col <- argv$sample_col
gene_id_type <- argv$id_type
fdr_cutoff <- argv$fdr
c1 <- argv$c1
c2 <- argv$c2
contrast <- c(cond_col, c1, c2)
n_cpus <- argv$cpus
covariate_formula <- argv$covariate_formula

register(MulticoreParam(workers = n_cpus))

# Reading the Annotation sample csv file
sampleAnno <- read_csv(argv$colData) |> filter(get(cond_col) %in% contrast[2:3])
# Reading the Count matrix tsv file
count_mat <- read_csv(argv$countData)

# if (gene_id_type == "ENSEMBL") {
#     count_mat <- count_mat |> mutate(gene_id = remove_ensg_version(gene_id))
# }

# ensg_to_genesymbol <- count_mat |> dplyr::select(gene_id, gene_name)
# ensg_to_desc <- AnnotationDbi::select(org.Hs.eg.db, count_mat$gene_id |> unique(), keytype = gene_id_type, columns = c("GENENAME")) |>
#     distinct(across(!!gene_id_type), .keep_all = TRUE)

count_mat <- count_mat |>
    dplyr::select(c(gene_id, sampleAnno[[sample_col]])) |>
    column_to_rownames("gene_id") |>
    round() # salmon does not necessarily contain integers

design_formula <- as.formula(paste0("~", cond_col, " ", covariate_formula))

# Using this subclass to store input values,intermediate calculations and result of DE analysis
dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = sampleAnno,
    design = design_formula
)

## keep only genes where we have >= 10 reads per samplecondition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep, ]

# save filtered count file
write_tsv(counts(dds) |> as_tibble(rownames = "gene_id"), file.path("./", paste0(prefix, "_detectedGenesRawCounts_min_10_reads_in_one_condition.tsv")))

# save normalized filtered count file
dds <- estimateSizeFactors(dds)
write_tsv(counts(dds, normalized = TRUE) |> as_tibble(rownames = "gene_id"), file.path("./", paste0(prefix, "_detectedGenesNormalizedCounts_min_10_reads_in_one_condition.tsv")))

# Set the reference to the contrast level 2 (baseline) given by the --c2 option
dds[[cond_col]] <- relevel(dds[[cond_col]], contrast[[3]])

# run DESeq
dds <- DESeq(dds, parallel = (n_cpus > 1))

### IHW
# use of IHW for p value adjustment of DESeq2 results
resIHW <- results(dds, filterFun = ihw, contrast = contrast) |>
    as_tibble(rownames = "gene_id") |>
    # left_join(ensg_to_genesymbol) |>
    # left_join(ensg_to_desc, by = c("gene_id" = gene_id_type)) |>
    # dplyr::rename(genes_description = GENENAME) |>
    arrange(pvalue)

write_tsv(resIHW |> as_tibble(rownames = "gene_id"), file.path("./", paste0(prefix, "_DESeq2_result.tsv")))
summary(resIHW)
ngenes_cut <- sum(resIHW$padj < fdr_cutoff, na.rm = TRUE)
print(paste0("Number of genes under the specified cutoff: ", ngenes_cut))
