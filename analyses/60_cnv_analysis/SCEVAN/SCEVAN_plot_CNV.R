# Load libraries
library(rstatix)
library(stats)
library(ggpubr)
library(stringr)
library(plotly)
library(tidyr)
library(dplyr)
library(circlize)
library(ComplexHeatmap)
library(svglite)

# Load functions
loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

kos <- function(x){
  # read data and set column names
  data <- loadRData(x)
  # get the sample name from the file path
  sampleName <- str_extract(x, "full_atlas_merged_.+_count_mtx_annot.RData")
  sampleName <- gsub("_count_mtx_annot.RData", "", sampleName)
  sampleName <- gsub("full_atlas_annotated_", "", sampleName)
  data$sample <- sampleName
  return(data)
}

in_path <- "./data/"

low_fdr <- read.csv(paste0(in_path,"cnv_lm_res_filtered_abs_diff_gt_0.01_fdr_lt_0.1_adjusted.csv"))

## Create the DF fo rall groups
circos_all_groups <- data.frame(
  chromosome = low_fdr$chromosome,
  start = low_fdr$start,
  end = low_fdr$end,
  ID_mean = low_fdr$segmean_desert,
  mixed_mean = low_fdr$segmean_B,
  M_mean = low_fdr$segmean_M,
  T_mean = low_fdr$segmean_T
)

circos_all_groups <- circos_all_groups %>%
  dplyr::arrange(start) %>% 
  dplyr::arrange(chromosome)

circos_all_groups <- unique(circos_all_groups)

### Gene annotation labels
# Get reference genomic positions from SCEVAN annotation files
groups <- c("Immune_desert", "mixed", "M", "T")

for(group in groups){
  meta_files <- Sys.glob(paste0(in_path,group,"/ann_mtx/*"))
  
  cnMetadata <- lapply(meta_files, kos)
  
  cnMetadata <- do.call("rbind", cnMetadata)
  
  assign(paste("metaData_", group, sep = ""), data.frame(
    "chromosome" = cnMetadata$seqnames, 
    "start" = cnMetadata$start,
    "end" = cnMetadata$end,
    "gene_name" = cnMetadata$gene_name
  ))
}

meta_merged <- rbind(metaData_Immune_desert, metaData_mixed, metaData_M, metaData_T)

final_metaData <- unique(meta_merged)

SCEVAN_df_annotated <- merge(circos_all_groups, final_metaData, by = c("chromosome", "start", "end"))

################################### CIRCOS ######################
# Build ref bed for ideogram
ref_bed <- read.cytoband(species = "hg38")

ref_df <- ref_bed$df
ref_df <- ref_df[ref_df$V1 != "chrX", ]
ref_df <- ref_df[ref_df$V1 != "chrY", ]

# Add "chr" to chromosomes - for plotting
circos_all_groups$chromosome <- paste("chr", circos_all_groups$chromosome , sep="")

# Build the circos DFs for each immune group
clean_circos_id <- data.frame(
  chr = circos_all_groups$chromosome,
  start = as.numeric(circos_all_groups$start),
  end = as.numeric(circos_all_groups$end),
  value1 = circos_all_groups$ID_mean
)


clean_circos_mixed <- data.frame(
  chr = circos_all_groups$chromosome,
  start = as.numeric(circos_all_groups$start),
  end = as.numeric(circos_all_groups$end),
  value1 = circos_all_groups$mixed_mean
)


clean_circos_m <- data.frame(
  chr = circos_all_groups$chromosome,
  start = as.numeric(circos_all_groups$start),
  end = as.numeric(circos_all_groups$end),
  value1 = circos_all_groups$M_mean
)

clean_circos_t <- data.frame(
  chr = circos_all_groups$chromosome,
  start = as.numeric(circos_all_groups$start),
  end = as.numeric(circos_all_groups$end),
  value1 = circos_all_groups$T_mean
)

### Circos
## Open file stream
## PNG
# png(file.path(paste0(in_path,"circos_SCEVAN_NSCLC.png")), 
#     width = 2800, height = 1600, units = "px", bg = "white",  res = 300)
# 
## PDF
# pdf(file.path(paste0(in_path,"circos_SCEVAN_NSCLC.pdf")), 
#     width = 9.3, height = 5.3, bg = "white",  pointsize = 10)

## SVG
svglite(file.path(paste0(in_path,"circos_SCEVAN_NSCLC.svg")), 
        width = 9.3, height = 5.3, bg = "white",  pointsize = 10)

# Start ideogram
circos.clear()
circos.par(gap.after=3)   #Increase gap size
circos.initializeWithIdeogram(ref_df, 
                              track.height = 0.005,
                              ideogram.height = 0.02)

# Loop through immune groups to create circos tracks
im_groups <- c("clean_circos_id", "clean_circos_mixed", "clean_circos_m", "clean_circos_t")

for (grp in im_groups) {
  circos.genomicTrack(get(grp), ylim = c(min(get(grp)[["value1"]]), max(get(grp)[["value1"]])), track.height = 0.09,
                      panel.fun = function(region, value, ...)  {
                        for(h in seq(min(get(grp)[["value1"]]), max(get(grp)[["value1"]]), by = 0.2)) {
                          circos.lines(CELL_META$cell.xlim, c(h, h), lty = 3, col = "#AAAAAA")
                        }
                        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 3, col = "#888888")
                        
                        circos.genomicPoints(region, value,
                                             col = ifelse(value[[1]] > 0, "#E41A1C", "#377EB8"),
                                             pch = 20, cex = 0.5)
                      })
}

### Add track labels
circos.text(sector.index="chr1",track.index = 2,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-2*max(get.cell.meta.data("cell.ylim")), labels = "A",facing = "clockwise", 
            niceFacing = TRUE, adj = c(3,0.1),cex=0.9)

circos.text(sector.index="chr1",track.index = 3,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-2*max(get.cell.meta.data("cell.ylim")), labels = "B",facing = "clockwise", 
            niceFacing = TRUE, adj = c(-2,0.3),cex=0.9)

circos.text(sector.index="chr1",track.index = 4,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-2*max(get.cell.meta.data("cell.ylim")), labels = "C",facing = "clockwise", 
            niceFacing = TRUE, adj = c(-1.5,0.3),cex=0.9)

circos.text(sector.index="chr1",track.index = 5,get.cell.meta.data("cell.xlim")-mean(get.cell.meta.data("cell.xlim"))/2,
            get.cell.meta.data("cell.ylim")-2*max(get.cell.meta.data("cell.ylim")), labels = "D",facing = "clockwise", 
            niceFacing = TRUE, adj = c(-3.5,0.4),cex=0.9)


# Use the genes of interest DF to add gene labels
shared_goi <- data.frame(
  gene_name = c("TMEM42", "ZDHHC3", "EXOSC7", "SACM1L", "KLHL18", "ELP6", "TCTA", "MAPKAPK3", "DCP1A", "ESF1", 
                "SNRPB2", "DSTN", "RRBP1", "SNX5", "MGME1", "SEC23B", "RIN2",
                "IGF2BP2", "SPON2", "PDCD10", "TNFSF10", "IL1RAP")
)

diff_genes <- merge(shared_cna, shared_goi, by = "gene_name")
diff_genes$chromosome <- paste("chr", diff_genes$chromosome , sep="")

diff_genes_final <- data.frame(
  chromosome = diff_genes$chromosome,
  start = as.numeric(diff_genes$start),
  end = as.numeric(diff_genes$end),
  gene_name = diff_genes$gene_name
)

circos.genomicLabels(diff_genes_final, labels.column = 4, side = "inside")

# Add legend
lgd_points = Legend(at = c("Gain", "Loss"), type = "points", 
                    legend_gp = gpar(col = c("#E41A1C", "#377EB8")), title_position = "topcenter", 
                    title = "CNA")

lgd_tracks = Legend(at = c("A: Immune Desert", "B: Mixed", "C: M infiltrated", "D: T infiltrated"), type = "grid", 
                    title_position = "topcenter", 
                    title = "Tracks")

lgd_list_vertical = packLegend(lgd_tracks, lgd_points)

circle_size = unit(1.45, "snpc")
draw(lgd_tracks, x = unit(1.39, "snpc"), y= unit(65, "mm"), just = "left")
draw(lgd_points, x = unit(1.41, "snpc"), y= unit(45, "mm"), just = "left")

circos.clear()
dev.off()