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
p <- add_argument(p, "--prefix", help = "Prefix of result file", default = "test")
p <- add_argument(p, "--sample_col", help = "Column in sample annotation that contains the sample names", default = "sample")
p <- add_argument(p, "--fdr", help = "False discovery rate for GO analysis and volcano plots", default = 0.1)
p <- add_argument(p, "--c1", help = "Contrast level 1 (perturbation)", default = "grpA")
p <- add_argument(p, "--c2", help = "Contrast level 2 (baseline)", default = "grpB")
p <- add_argument(p, "--cpus", help = "Number of cpus", default = 8)
p = add_argument(p, "--sum2zero", help = "Perform an all-vs-all comparison with sum2zero coding. If this flag is enabled the `--c1` and `--c2` flags are ignored.", flag=TRUE)

# Parse the command line arguments
argv <- parse_args(p)

# Function for removing the ENSG version
remove_ensg_version <- function(x) gsub("\\.[0-9]*$", "", x)

# # For testing only
# argv = environment()
# argv$results_dir = "/home/sturm/Downloads/sum2zero/"
# argv$prefix = "test"
# argv$cond_col = "condition"
# argv$sample_col = "sample"
# argv$fdr = 0.1
# argv$sum2zero = TRUE
# argv$c1 = NULL
# argv$c2 = NULL
# argv$cpus = 8
# argv$covariate_formula = ""
# argv$colData = "/data/projects/2020/Pircher-scRNAseq-lung/30_downstream_analyses/de_analysis/luad_lusc/pseudobulk_by_cell_type/adata_primary_tumor_tumor_cells_samplesheet.csv"
# argv$countData = "/data/projects/2020/Pircher-scRNAseq-lung/30_downstream_analyses/de_analysis/luad_lusc/pseudobulk_by_cell_type/adata_primary_tumor_tumor_cells_counts.csv"


# # For testing only
# argv = environment()
# argv$results_dir = "/home/sturm/Downloads/sum2zero/"
# argv$prefix = "test"
# argv$cond_col = "immune_infiltration"
# argv$sample_col = "sample"
# argv$fdr = 0.1
# argv$sum2zero = TRUE
# argv$c1 = NULL
# argv$c2 = NULL
# argv$cpus = 8
# argv$covariate_formula = "+ dataset"
# argv$colData = "/data/projects/2020/Pircher-scRNAseq-lung/30_downstream_analyses/de_analysis/immune_infiltration/pseudobulk_by_cell_type/adata_primary_tumor_tumor_cells_samplesheet.csv"
# argv$countData = "/data/projects/2020/Pircher-scRNAseq-lung/30_downstream_analyses/de_analysis/immune_infiltration/pseudobulk_by_cell_type/adata_primary_tumor_tumor_cells_counts.csv"


results_dir <- argv$resDir
prefix <- argv$prefix
cond_col <- argv$cond_col
sample_col <- argv$sample_col
fdr_cutoff <- argv$fdr
sum2zero = argv$sum2zero
if(!sum2zero) {
  c1 <- argv$c1
  c2 <- argv$c2
} else {
  c1 = NULL
  c2 = NULL
}
n_cpus <- argv$cpus
covariate_formula <- argv$covariate_formula

register(MulticoreParam(workers = n_cpus))

# Reading the Annotation sample csv file
sampleAnno <- read_csv(argv$colData)
if(!sum2zero) {
  sampleAnno = sampleAnno |> filter(get(cond_col) %in% c(c1, c2))
}
# Reading the Count matrix tsv file
count_mat <- read_csv(argv$countData)

count_mat <- count_mat |>
    dplyr::select(c(gene_id, sampleAnno[[sample_col]])) |>
    column_to_rownames("gene_id") |>
    round() # salmon does not necessarily contain integers

# mock data for testing only! 
# sampleAnno = tibble(sample=c("s1", "s2", "s3", "s4", "s5", "s6"), condition=c("A", "A", "B", "B", "C", "C"))
# count_mat = t(as.matrix(data.frame(gene1=c(4,6, 9,10, 1000,20000), gene2=c(4,6, 9,10, 2,1))))
# colnames(count_mat) = sampleAnno$sample

design_formula <- as.formula(paste0("~", cond_col, " ", covariate_formula))

if(sum2zero) {
  design_mat = model.matrix(design_formula, contrasts.arg = structure(as.list("contr.sum"), names=cond_col), data=sampleAnno)  
} else {
  design_mat = design_formula
}

# Using this subclass to store input values,intermediate calculations and result of DE analysis
dds <- DESeqDataSetFromMatrix(
    countData = count_mat,
    colData = sampleAnno,
    design = design_mat 
)

## keep only genes where we have >= 10 reads per samplecondition in total
keep <- rowSums(counts(collapseReplicates(dds, dds[[cond_col]]))) >= 10
dds <- dds[keep, ]

# save filtered count file
write_tsv(counts(dds) |> as_tibble(rownames = "gene_id"), file.path("./", paste0(prefix, "_detectedGenesRawCounts_min_10_reads_in_one_condition.tsv")))

# save normalized filtered count file
dds <- estimateSizeFactors(dds)
write_tsv(counts(dds, normalized = TRUE) |> as_tibble(rownames = "gene_id"), file.path("./", paste0(prefix, "_detectedGenesNormalizedCounts_min_10_reads_in_one_condition.tsv")))

# run DESeq
dds <- DESeq(dds, parallel = (n_cpus > 1))

if(sum2zero) {
  # order needs to be the one of the levels of the factor (same as for contrast matrix)
  unique_conditions = levels(as.factor(sampleAnno[[cond_col]]))
  n_unique = length(unique_conditions)
  # with sum2zero we test that a single coefficient != 0
  # a coefficient corresponds to the difference from the overall mean
  # the intercept correponds to the overall mean
  contr_mat = diag(n_unique - 1) 
  # make list with one contrast per item
  contrasts = lapply(seq_len(n_unique - 1), function(i) { contr_mat[, i] }) 
  # the above added n-1 comparisons, we need to construct the last (redundant) one manually
  contrasts = append(contrasts, list(-apply(contr_mat, MARGIN = 1, sum) / (n_unique - 1)))
  # pad end of vector with zeros (required if there are covariates in the design).
  # we can assume that the "condition columns" always come at the front since
  # it is the first argument of the formula
  contrasts = lapply(contrasts, function(x) {
    c(0, x, rep.int(0, length(resultsNames(dds)) - n_unique))
  })
  # set names of contrasts
  names(contrasts) = unique_conditions
} else {
  contrasts = list(c(cond_col, c1, c2))
  names(contrasts) = sprintf("%s_vs_%s", c1, c2)
}

### IHW
# use of IHW for p value adjustment of DESeq2 results
resIHW = lapply(names(contrasts), function(name) {
  contrast = contrasts[[name]]
  results(dds, filterFun = ihw, contrast = contrast) |>
    as_tibble(rownames = "gene_id") |>
    mutate(comparison = name) |>
    arrange(pvalue)
}) |> bind_rows()


write_tsv(resIHW, file.path("./", paste0(prefix, "_DESeq2_result.tsv")))
summary(resIHW)
ngenes_cut <- sum(resIHW$padj < fdr_cutoff, na.rm = TRUE)
print(paste0("Number of genes under the specified cutoff: ", ngenes_cut))
