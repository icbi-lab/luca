library("tibble")
library("readr")
library("dplyr")
library("tidyr")
library("ComplexHeatmap")
library("circlize")
library("GenomicDataCommons")
library("foreach")
library("doParallel")
library("clusterProfiler")
library("org.Hs.eg.db")
library("biomaRt")
library("BiocParallel")
library("ggplot2")
library("reticulate")
library("svglite")

library("conflicted")

conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("paste", "base")
conflict_prefer("rename", "dplyr")

## register cores for parallel computing
registerDoParallel(cores=32)
register(MulticoreParam(workers=32))

## result data and python script dirs
data_dir = "./data"
bin_dir = "./bin"

## circos plot output format
circos_pdf = TRUE
circos_svg = FALSE

## celltype fraction files
TCGA_LUAD_ctf <- file.path(data_dir, "TCGA_LUAD_quanTIseq_lsei_TIL10_cellfractions.txt")
TCGA_LUSC_ctf <- file.path(data_dir, "TCGA_LUSC_quanTIseq_lsei_TIL10_cellfractions.txt")

## set python conda env
options(reticulate.conda_binary = "/usr/local/bioinf/conda/condabin/conda")
use_condaenv("CNAutils")

# get python path
python_bin = conda_python("CNAutils")
python_path =gsub("\\/python$", "",python_bin, perl = TRUE)
env_path = paste0(python_path,":",Sys.getenv("PATH"))
Sys.setenv(PATH = env_path)

## read quanTIseq cell type fraction data from TCIA
cell_types_LUAD <- read_tsv(TCGA_LUAD_ctf) |> column_to_rownames(var="CellType")
cell_types_LUSC <- read_tsv(TCGA_LUSC_ctf) |> column_to_rownames(var="CellType")

# join LUAD and LUSC and convert to matrix
cell_types_NSCLC <- as.matrix(bind_cols(cell_types_LUAD, cell_types_LUSC))

## number of clusters to create
nr_cl <- 4

## calculate patient-wise z-score over celltypes
cell_types_NSCLC_z <- t(apply(cell_types_NSCLC, 1, scale))
colnames(cell_types_NSCLC_z) <- colnames(cell_types_NSCLC)


## colorbar settings for heatmap
max_cb <- min(abs(min(cell_types_NSCLC_z)), abs(max(cell_types_NSCLC_z)))
col_fun = colorRamp2(c(-max_cb, 0, max_cb), c("blue", "white", "red"))

## cluster celltype fractions (z-scores)
# set seed
set.seed(13)
HM <- Heatmap(cell_types_NSCLC_z,
        show_column_names = FALSE,
        cluster_rows = FALSE,
        column_km = nr_cl,
        column_km_repeats = 100,
        column_dend_height = unit(30, "mm"),
        col = col_fun,
        heatmap_legend_param = list(title = "cell type fraction z-score")
)

png(file.path(data_dir, "celltype_fractions_clustered.png"), width = 2200, height = 1000, units = "px", pointsize = 6,
    bg = "white",  res = 300)
HM <- ComplexHeatmap::draw(HM)
dev.off()


## extract clusters
co <- column_order(HM)
# cluster 1 is immune desert ID
immune_desert <- colnames(cell_types_NSCLC_z)[co[["1"]]]

# clusters 2-4 are immune infiltrated
immune_infiltrated <- colnames(cell_types_NSCLC_z)[c(co[["2"]], co[["3"]], co[["4"]])]

# to_export = c(rep_len("desert", length(immune_desert)), rep_len("infiltrated", length(immune_infiltrated)))
# names(to_export) = c(immune_desert, immune_infiltrated)
# to_export_df = tibble(patient_strat = to_export, TCGA_patient_barcode=names(patient_strat))
# write_tsv(to_export_df, "../../../tables/tcga/quantiseq_patient_stratification_didi.tsv")

## Download CNA data from GDC
# segments
cna_id_manifest = files() %>%
  GenomicDataCommons::filter(cases.submitter_id == immune_desert) %>%
  GenomicDataCommons::filter(cases.samples.sample_type == 'primary tumor') %>%
  GenomicDataCommons::filter(data_type == 'Masked Copy Number Segment') %>%
  manifest()

cna_if_manifest = files() %>%
  GenomicDataCommons::filter(cases.submitter_id == immune_infiltrated) %>%
  GenomicDataCommons::filter(cases.samples.sample_type == 'primary tumor') %>%
  GenomicDataCommons::filter(data_type == 'Masked Copy Number Segment') %>%
  manifest()

# genes
cna_genes_id_manifest = files() %>%
  GenomicDataCommons::filter(cases.submitter_id == immune_desert) %>%
  GenomicDataCommons::filter(cases.samples.sample_type == 'primary tumor') %>%
  GenomicDataCommons::filter(data_type == 'Gene Level Copy Number') %>%
  manifest()

cna_genes_if_manifest = files() %>%
  GenomicDataCommons::filter(cases.submitter_id == immune_infiltrated) %>%
  GenomicDataCommons::filter(cases.samples.sample_type == 'primary tumor') %>%
  GenomicDataCommons::filter(data_type == 'Gene Level Copy Number') %>%
  manifest()



## get filenames
cna_id_fnames = gdcdata(cna_id_manifest$id,progress=FALSE)
cna_if_fnames = gdcdata(cna_if_manifest$id,progress=FALSE)
cna_genes_id_fnames = gdcdata(cna_genes_id_manifest$id,progress=FALSE)
cna_genes_if_fnames = gdcdata(cna_genes_if_manifest$id,progress=FALSE)

## save GDC data from cache dir into data dir
if(! file.exists(data_dir)) {
  dir.create(data_dir)
}

data_dir_id = file.path(data_dir, "immune_desert")
if(! file.exists(data_dir_id)) {
  dir.create(data_dir_id)
}

data_dir_if = file.path(data_dir, "immune_infiltrated")
if(! file.exists(data_dir_if)) {
  dir.create(data_dir_if)
}

data_dir_id_genes = file.path(data_dir, "immune_desert_genes")
if(! file.exists(data_dir_id_genes)) {
  dir.create(data_dir_id_genes)
}

data_dir_if_genes = file.path(data_dir, "immune_infiltrated_genes")
if(! file.exists(data_dir_if_genes)) {
  dir.create(data_dir_if_genes)
}

for (d in names(cna_id_fnames)) {
  file.copy(list.files(file.path("~/.cache/GenomicDataCommons/", d), full.names = TRUE), data_dir_id, recursive = TRUE) 
}

for (d in names(cna_if_fnames)) {
  file.copy(list.files(file.path("~/.cache/GenomicDataCommons/", d), full.names = TRUE), data_dir_if, recursive = TRUE) 
}

for (d in names(cna_genes_id_fnames)) {
  file.copy(list.files(file.path("~/.cache/GenomicDataCommons/", d), full.names = TRUE), data_dir_id_genes, recursive = TRUE) 
}

for (d in names(cna_genes_if_fnames)) {
  file.copy(list.files(file.path("~/.cache/GenomicDataCommons/", d), full.names = TRUE), data_dir_if_genes, recursive = TRUE) 
}


## map GCD file ids to  TCGA case id
# Function to map GDC file ids to TCGA case ids (case.submitter_id)

TCGAtranslateID = function(file_ids, legacy = FALSE) {
    info = files(legacy = legacy) %>%
        GenomicDataCommons::filter( ~ file_id %in% file_ids) %>%
        GenomicDataCommons::select('cases.submitter_id') %>%
        results_all()
    # The mess of code below is to extract TCGA barcodes
    # id_list will contain a list (one item for each file_id)
    # of TCGA barcodes of the form 'TCGA-XX-YYYY-ZZZ'
    id_list = lapply(info$cases,function(a) {
        a[[1]][[1]][[1]]})
    # so we can later expand to a data.frame of the right size
    barcodes_per_file = sapply(id_list,length)
    # And build the data.frame
    return(data.frame(file_id = rep(ids(info),barcodes_per_file),
                      submitter_id = unlist(id_list)))
}


## do mapping and return distinct submitter_id (multiple files per patients may be present)
cna_id_fid_bc <- TCGAtranslateID(names(cna_id_fnames)) |> distinct_at("submitter_id")
cna_if_fid_bc <- TCGAtranslateID(names(cna_if_fnames)) |> distinct_at("submitter_id")
cna_genes_id_fid_bc <- TCGAtranslateID(names(cna_genes_id_fnames)) |> distinct_at("submitter_id")
cna_genes_if_fid_bc <- TCGAtranslateID(names(cna_genes_if_fnames)) |> distinct_at("submitter_id")

## Make and export filename TCGA barcode map for immune desert (id): genes
cna_genes_id_bc_fnames <- cna_genes_id_fnames[rownames(cna_genes_id_fid_bc)] %>%
  bind_cols(cna_genes_id_fid_bc$submitter_id) %>%
  as_tibble() %>%
  set_colnames(c("fname", "bc"))

## Make and export filename TCGA barcode map for immune infiltrated (if): genes
cna_genes_if_bc_fnames <- cna_genes_if_fnames[rownames(cna_genes_if_fid_bc)] %>%
  bind_cols(cna_genes_if_fid_bc$submitter_id) %>%
  as_tibble() %>%
  set_colnames(c("fname", "bc"))


## genes
write_tsv(cna_genes_id_bc_fnames, file.path(data_dir, "cna_genes_id_files.txt"), col_names = FALSE)
write_tsv(cna_genes_if_bc_fnames, file.path(data_dir, "cna_genes_if_files.txt"), col_names = FALSE)


## Make and export filename TCGA barcode map for immune desert (id): segments
cna_id_bc_fnames <- cna_id_fnames[rownames(cna_id_fid_bc)] %>%
  bind_cols(cna_id_fid_bc$submitter_id) %>%
  as_data_frame() %>%
  set_colnames(c("fname", "bc"))

## Make and export filename TCGA barcode map for immune infiltrated (if): segments
cna_if_bc_fnames <- cna_if_fnames[rownames(cna_if_fid_bc)] %>%
  bind_cols(cna_if_fid_bc$submitter_id) %>%
  as_data_frame() %>%
  set_colnames(c("fname", "bc"))

## segments to genes filename mapping
cna_seg_genes_id_bc_fnames <- cna_id_bc_fnames %>% inner_join(cna_genes_id_bc_fnames, by = "bc", suffix = c(".seg", ".genes"))
cna_seg_genes_if_bc_fnames <- cna_if_bc_fnames %>% inner_join(cna_genes_if_bc_fnames, by = "bc", suffix = c(".seg", ".genes"))

write_tsv(cna_seg_genes_id_bc_fnames, file.path(data_dir, "cna_seg_genes_id_files.txt"), col_names = FALSE)
write_tsv(cna_seg_genes_if_bc_fnames, file.path(data_dir, "cna_seg_genes_if_files.txt"), col_names = FALSE)

## make new bc fname map with new gene level cna from seg files
# id
cna_genes_id_bc_fnames <- data.frame(fname=character(), bc=character())
sys_cmds <- data.frame(cmd=character())

id_genes_seg_dir = file.path(data_dir, "immune_desert_genes_seg")
if(! file.exists(id_genes_seg_dir)) {
  dir.create(id_genes_seg_dir)
}

for (i in 1:nrow(cna_seg_genes_id_bc_fnames)) {
  print(cna_seg_genes_id_bc_fnames$bc[i])
  outfile <- file.path(id_genes_seg_dir, paste0(cna_seg_genes_id_bc_fnames$bc[i], "_gene_seg.tsv"))
  sys_cmd <- paste(python_bin,
                 file.path(bin_dir, "make_intersected_bed.py"),
                 "--gene_cnv", cna_seg_genes_id_bc_fnames$fname.genes[i],
                 "--seg_cnv", cna_seg_genes_id_bc_fnames$fname.seg[i],
                 "--outfile", outfile,
                 sep = " "
                 )
  sys_cmds[nrow(sys_cmds)+1,] <- sys_cmd

  cna_genes_id_bc_fnames[nrow(cna_genes_id_bc_fnames)+1, ] <- c(outfile, cna_seg_genes_id_bc_fnames$bc[i])
}


# secure this, will run for long time
if (! exists(".cna_genes_id_done") & ! file.exists(file.path(data_dir, ".cna_genes_id_done"))) {
  foreach(i=1:nrow(sys_cmds), .combine=rbind) %dopar% { system(sys_cmds[i, ])}
  write_tsv(cna_genes_id_bc_fnames, file.path(data_dir, "cna_genes_id_files.txt"), col_names = FALSE)
  .cna_genes_id_done <- 1
  file.create(".cna_genes_id_done")
}


## run python script that pastes all "copy_number" columns for id into a single table and calculates some stats
sys_cmd <- paste(python_bin,
                 file.path(bin_dir, "paste_cn.py"),
                 "6",
                 file.path(data_dir, "cna_genes_id_files.txt"),
                 "0.1",
                 ">",
                 file.path(data_dir, "cna_genes_id.tsv"), sep = " ")

system(sys_cmd)


# if
cna_genes_if_bc_fnames <- data.frame(fname=character(), bc=character())
sys_cmds <- data.frame(cmd=character())

if_genes_seg_dir = file.path(data_dir, "immune_infiltrated_genes_seg")
if(! file.exists(if_genes_seg_dir)) {
  dir.create(if_genes_seg_dir)
}

for (i in 1:nrow(cna_seg_genes_if_bc_fnames)) {
  print(cna_seg_genes_if_bc_fnames$bc[i])
  outfile <- file.path(if_genes_seg_dir, paste0(cna_seg_genes_if_bc_fnames$bc[i], "_gene_seg.tsv"))
  sys_cmd <-  paste(python_bin,
                   file.path(bin_dir, "make_intersected_bed.py"),
                   "--gene_cnv", cna_seg_genes_if_bc_fnames$fname.genes[i],
                   "--seg_cnv", cna_seg_genes_if_bc_fnames$fname.seg[i],
                   "--outfile", outfile,
                   sep = " "
  )
  sys_cmds[nrow(sys_cmds)+1,] <- sys_cmd

  cna_genes_if_bc_fnames[nrow(cna_genes_if_bc_fnames)+1, ] <- c(outfile, cna_seg_genes_if_bc_fnames$bc[i])
}

# secure this, will run for long time
if (! exists(".cna_genes_if_done") & ! file.exists(file.path(data_dir, ".cna_genes_if_done"))) {
  foreach(i=1:nrow(sys_cmds), .combine=rbind) %dopar% { system(sys_cmds[i, ])}
  write_tsv(cna_genes_if_bc_fnames, file.path(data_dir, "cna_genes_if_files.txt"), col_names = FALSE)
  .cna_genes_if_done <- 1
  file.create(".cna_genes_if_done")
}

## run python script that pastes all "copy_number" columns for if into a single table and calculates some stats
sys_cmd <- paste(python_bin,
                 file.path(bin_dir, "paste_cn.py"),
                 "6",
                 file.path(data_dir, "cna_genes_if_files.txt"),
                 "0.1",
                 ">",
                 file.path(data_dir, "cna_genes_if.tsv"), sep = " ")

system(sys_cmd)


#
## wilcox.test single genes to get significant diff ID vs IF
#
cna_genes_id_all <- read_tsv(file.path(data_dir, "cna_genes_id.tsv")) %>%
  select(c(gene_name, gene_id, starts_with("TCGA")))
  
cna_genes_if_all <- read_tsv(file.path(data_dir, "cna_genes_if.tsv")) %>%
  select(c(gene_name, gene_id, starts_with("TCGA")))

# make matrix: calculates faster
cna_genes_id_all_m <- cna_genes_id_all[,3:ncol(cna_genes_id_all)] %>% as.matrix()
cna_genes_if_all_m <- cna_genes_if_all[,3:ncol(cna_genes_if_all)] %>% as.matrix()

# data frame for stat
cna_genes_stat <- data.frame("gene_id" = character(), "gene_name" = character(), "p.val" = numeric(), "p.adj" = numeric())

## get nr of patients
id_n <- length(rownames(cna_genes_id_fid_bc))
if_n <- length(rownames(cna_genes_if_fid_bc))

# run wilcox.tests
for (i in 1:nrow(cna_genes_id_all)) {
  cna_id <- cna_genes_id_all_m[i,]
  cna_if <- cna_genes_if_all_m[i,]
  
  if (sum(is.na(cna_id)) != id_n & sum(is.na(cna_if)) != if_n) {
    # we do not use t.test, data not always normal distributed
    # cna_gene_stats <- t.test(cna_id, cna_if, na.action="na.exclude")
    # use non-parametric test
    cna_gene_stats <- wilcox.test(cna_id, cna_if, na.action="na.exclude")
    pval <- cna_gene_stats$p.value
    print(i)
  } else {
    pval <- NaN
  }
  cna_genes_stat[i,] <- c(cna_genes_id_all$gene_id[i], cna_genes_id_all$gene_name[i], pval, NaN)
}

## get biotype of genes, to filter for protein coding ones

# init ensembl
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

# function to remove versions
remove_ensg_version = function(x) gsub("\\.[0-9]*$", "", x)

# genes to lookup in ensembl
goi <- remove_ensg_version(cna_genes_stat$gene_id)

# use biomart to get biotype
goids = getBM(attributes = c('ensembl_gene_id', 'gene_biotype'), 
              filters = 'ensembl_gene_id', 
              values = goi, 
              mart = ensembl)

## filter for protein coding genes
cna_genes_stat <- cna_genes_stat %>%
  mutate(gene_id = remove_ensg_version(gene_id)) %>%
  inner_join(goids, by = c("gene_id" = "ensembl_gene_id")) %>%
  filter(gene_biotype == "protein_coding")

# adjust p.value
cna_genes_stat$p.adj <- p.adjust(cna_genes_stat$p.val, "BH")


#### Start circos wilcox.test

## genes, calculate means and filter
# id
cna_genes_id_filtered_pc <- read_tsv(file.path(data_dir, "cna_genes_id.tsv")) %>%
  select(gene_id, gene_name, chromosome, start, end, starts_with("TCGA")) %>%
  mutate(gene_id = remove_ensg_version(gene_id)) %>%
  gather(SampleID, mean_cn, starts_with("TCGA")) %>%
  group_by(gene_id, gene_name, chromosome, start) %>%
  select(-SampleID) %>%
  summarise_each(funs(mean)) %>%
  inner_join(goids, by = c("gene_id" = "ensembl_gene_id")) %>%
  filter(gene_biotype == "protein_coding")


# if
cna_genes_if_filtered_pc <- read_tsv(file.path(data_dir, "cna_genes_if.tsv")) %>%
  select(gene_id, gene_name, chromosome, start, end, starts_with("TCGA")) %>%
  mutate(gene_id = remove_ensg_version(gene_id)) %>%
  gather(SampleID, mean_cn, starts_with("TCGA")) %>%
  group_by(gene_id, gene_name, chromosome, start) %>%
  select(-SampleID) %>%
  summarise_each(funs(mean)) %>%
  inner_join(goids, by = c("gene_id" = "ensembl_gene_id")) %>%
  filter(gene_biotype == "protein_coding")


# join id + if + and p-vals and filter
cna_genes_filtered_pc <- cna_genes_id_filtered_pc %>% 
  full_join(cna_genes_if_filtered_pc, by = c("gene_id", "gene_name", "chromosome", "start", "end", "gene_biotype"), suffix = c(".id", ".if")) %>%
  full_join(cna_genes_stat, by = c("gene_id", "gene_name", "gene_biotype")) %>%
  filter(
    ((abs(mean_cn.id) > 0.1 | abs(mean_cn.if > 0.1)))
    & (p.adj < 1e-4)
  )
dim(cna_genes_filtered_pc)[1]


## Create circos

# Build ref bed for ideogram without chrY and chrY
ref_bed <- read.cytoband(species = "hg38")
ref_df <- ref_bed$df |> filter(V1 != "chrX", V1 != "chrY")

# id only plot genes with abs(mean CNA) > 0.1 otherwise set 0.0 (= not CNA)
value1 <- ifelse((is.na(cna_genes_filtered_pc$mean_cn.id) | abs(cna_genes_filtered_pc$mean_cn.id) < 0.1), 0, cna_genes_filtered_pc$mean_cn.id)
clean_circos_tt_id <- data.frame(
  chr = cna_genes_filtered_pc$chromosome,
  start = cna_genes_filtered_pc$start,
  end = cna_genes_filtered_pc$end,
  value1 = value1
)

# if only plot genes with abs(mean CNA) > 0.1 otherwise set 0.0 (= not CNA)
value1 <- ifelse((is.na(cna_genes_filtered_pc$mean_cn.if) | abs(cna_genes_filtered_pc$mean_cn.if) < 0.1), 0, cna_genes_filtered_pc$mean_cn.if)
clean_circos_tt_if <- data.frame(
  chr = cna_genes_filtered_pc$chromosome,
  start = cna_genes_filtered_pc$start,
  end = cna_genes_filtered_pc$end,
  value1 = value1
)

## plot circos and save to png
if (circos_pdf == TRUE) {
    pdf(file.path(data_dir, "circos_TCGA_NSCLC_IDvsIF_3.pdf"), width = 9.3, height = 5.3, bg = "white",  pointsize = 10)
} else if (circos_svg == TRUE) {
    svglite(file.path(data_dir, "circos_TCGA_NSCLC_IDvsIF_3.svg"), width = 9.3, height = 5.3, bg = "white",  pointsize = 10)
} else {
    png(file.path(data_dir, "circos_TCGA_NSCLC_IDvsIF.png"), width = 2800, height = 1600, units = "px", bg = "white",  res = 300)
}
circos.clear()

# Increase gap size
circos.par(gap.after=3)

# plot ideogram
circos.initializeWithIdeogram(ref_df, 
                              track.height = 0.005,
                              ideogram.height = 0.02)

# plot id
circos.genomicTrack(clean_circos_tt_id, ylim = c(-0.4, 0.4), track.height = 0.09,
                    panel.fun = function(region, value, ...)  {
                      for(h in seq(-0.4, 0.4, by = 0.2)) {
                        circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#AAAAAA")
                      }
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 3, col = "#888888")

                      circos.genomicPoints(region, value, 
                                           col = ifelse(value[[1]] > 0.1, "#E41A1C", ifelse(value[[1]] < -0.1, "#377EB8", "#888888")), 
                                           # col = ifelse(value[[1]] > 0, "#E41A1C", "#377EB8"), 
                                           pch = 20, cex = 0.3)
                    })

# plot if
circos.genomicTrack(clean_circos_tt_if, ylim = c(-0.4, 0.4), track.height = 0.09,
                    panel.fun = function(region, value, ...)  {
                      for(h in seq(-0.4, 0.4, by = 0.2)) {
                        circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#AAAAAA")
                      }
                      circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 3, col = "#888888")

                      circos.genomicPoints(region, value, 
                                           col = ifelse(value[[1]] > 0.1, "#E41A1C", ifelse(value[[1]] < -0.1, "#377EB8", "#888888")), 
                                           # col = ifelse(value[[1]] > 0, "#E41A1C", "#377EB8"), 
                                           pch = 20, cex = 0.3)
                    })

# Genes of interest to label in circos (may be changed)
gene_list_to_mark <- as_tibble(
  unique(c("TMEM42", "ZDHHC3", "EXOSC7", "SACM1L", "KLHL18", "ELP6", "TCTA", "MAPKAPK3", "DCP1A", "ESF1", "SNRPB2", "DSTN", "RRBP1", "SNX5", "MGME1", "SEC23B", "RIN2",
           "IGF2BP2", "SPON2", "PDCD10", "TNFSF10", "IL1RAP"))
  ) %>%
  rename("gene_name" = value) %>%
  inner_join(cna_genes_filtered_pc, by = c("gene_name")) %>%
  select(c(chromosome, start, end, gene_name))

# add genes of interest labels
circos.genomicLabels(gene_list_to_mark, labels.column = 4, side = "inside")

# add track labels
circos.text(sector.index="chr1",track.index = 2,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-2*max(get.cell.meta.data("cell.ylim")), labels = "A",facing = "clockwise", 
            niceFacing = TRUE, adj = c(1.5,0.3),cex=0.4)

circos.text(sector.index="chr1",track.index = 3,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-2*max(get.cell.meta.data("cell.ylim")), labels = "B",facing = "clockwise", 
            niceFacing = TRUE, adj = c(-1.5,0.3),cex=0.4)

# Add legend
lgd_points = Legend(at = c("Gain", "Loss"), type = "points", 
                    legend_gp = gpar(col = c("#E41A1C", "#377EB8")), title_position = "topcenter", 
                    title = "CNA")

lgd_tracks = Legend(at = c("A: Immune desert", "B: Immune infiltrated"), type = "grid", 
                    title_position = "topcenter", 
                    title = "Tracks")

lgd_list_vertical = packLegend(lgd_tracks, lgd_points)
draw(lgd_tracks, x = unit(1.39, "snpc"), y= unit(65, "mm"), just = "left")
draw(lgd_points, x = unit(1.41, "snpc"), y= unit(45, "mm"), just = "left")

circos.clear()
dev.off()

#### End circos wilcox.test

# save ID/IF t.test sig genes
write_tsv(cna_genes_filtered_pc, file.path(data_dir, "ID_IF_CNA_diff_signifcant_protein_coding.tsv"))


## ORA for t.test sig genes with abs(mean CNA) > 0.3

# ORA output settings
prefix = "ID_IF_diff_sig"
ora_dir = file.path(data_dir, "/ORA")
if(! file.exists(ora_dir)) {
  dir.create(ora_dir)
}

# filter genes with |mean_cn| > 0.3
cna_genes_filtered_pc_ora <- cna_genes_filtered_pc %>% 
  filter(
    ((abs(mean_cn.id) > 0.3 | abs(mean_cn.if > 0.3)))
  )
dim(cna_genes_filtered_pc_ora)[1]

# create HGNC to ENTREZ id map
hgnc_to_entrez = AnnotationDbi::select(
  org.Hs.eg.db, cna_genes_stat %>% 
  pull("gene_name") %>%
  unique(), keytype="SYMBOL", columns=c("ENTREZID")
  )

# make background gene list
universe = hgnc_to_entrez %>% pull("ENTREZID") %>% unique()

# map HGNC symbols to ENTREZ ids in filtered gene list for ORA
genes = cna_genes_filtered_pc_ora %>% 
  inner_join(hgnc_to_entrez, by=c("gene_name"="SYMBOL"))

# calculate log2(CN ratio) ID/IF
cna_val <- log2(2^genes$mean_cn.id / 2^genes$mean_cn.if)
names(cna_val) <- genes$ENTREZID

# function for plot saving
save_plot <- function(filename, p, width=NULL, height=NULL) {
  if (!is.null(width) && !is.null(height)) {
    ggsave(file.path(paste0(filename, ".png")), plot = p, width = width, height = height, bg = "white")
    ggsave(file.path(paste0(filename, ".svg")), plot = p, width = width, height = height)
  } else {
    ggsave(file.path(paste0(filename, ".png")), plot = p, bg = "white")
    ggsave(file.path(paste0(filename, ".svg")), plot = p)
  }
}

# Warmup GO database - work around https://github.com/YuLab-SMU/clusterProfiler/issues/207
._ = enrichGO(universe[1], OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", universe = universe)

# function to calculate heatmap dimensions for heatplots
get_heatplot_dims <- function(p) {
  nr_gene <- length(unique(p$data$Gene))
  nr_cat <- length(unique(p$data$categoryID))
  
  hp_width = min(nr_gene * 0.25, 40)
  hp_height = min(nr_cat * 0.25, 40)
  
  return(c(hp_width, hp_height))
}

# define ORA tests to be performed: KEGG, WikiPathways, GO
ora_tests = list(
  "KEGG" = function(genes, universe) {
    enrichKEGG(
      gene         = genes,
      universe     = universe,
      organism     = 'hsa',
      pvalueCutoff = 0.05
    )
  },
  "WikiPathway" = function(genes, universe) {
    enrichWP(
      gene = genes,
      universe     = universe,
      organism     = 'Homo sapiens',
      pvalueCutoff = 0.05
    )
  },
  "GO_BP" = function(genes, universe) {
    enrichGO(
      gene = genes,
      universe = universe,
      keyType = "ENTREZID",
      OrgDb = org.Hs.eg.db,
      ont = "BP",
      pAdjustMethod = "BH",
      qvalueCutoff = 0.05,
      minGSSize = 10
    )
  }
)

# run ORA tests
bplapply(names(ora_tests), function(ora_name) {
  message(paste0("Performing ", ora_name, " ORA-test..."))
  
  test_fun = ora_tests[[ora_name]]
  ora_res = test_fun(genes$ENTREZID, universe)
  ora_res = setReadable(ora_res, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  res_tab = as_tibble(ora_res@result)
  
  write_tsv(res_tab, file.path(ora_dir, paste0(prefix, "_ORA_", ora_name, ".tsv")))
  
  if (min(res_tab$p.adjust) < 0.05) {
    p = dotplot(ora_res, showCategory=40)
    
    save_plot(file.path(ora_dir, paste0(prefix, "_ORA_", ora_name, "_dotplot")), p, width = 15, height = 10)
    
    p <- cnetplot(ora_res,
                  categorySize="pvalue",
                  showCategory = 5,
                  foldChange=cna_val,
                  vertex.label.font=6)
    
    save_plot(file.path(ora_dir, paste0(prefix, "_ORA_", ora_name, "_cnetplot")), p, width = 15, height = 12)
    
    p <- heatplot(ora_res, foldChange=cna_val, showCategory=40) +
      scale_fill_gradient2(midpoint=0, low="blue4", mid="white", high="red4" )
    hp_dims <- get_heatplot_dims(p)
    
    save_plot(file.path(ora_dir, paste0(prefix, "_ORA_", ora_name, "_heatplot")), p, width = hp_dims[1], height = hp_dims[2])
    
  } else {
    message(paste0("Warning: No significant enrichment in ", ora_name, " ORA analysis. "))
  }
})
