# BiocManager::install("cmapR")

library(cmapR)
library(dplyr)

compoundinfo <- read.delim(file = "compoundinfo_beta.txt")

# load metadata files
## -----------------------------------------------------------------------------

pertinfo <- read.delim(file = "LDS-1613/Data/pertinfo.txt")
colnames(pertinfo)
# [1] "pert_id"     "pert_mfc_id" "cmap_name"   "pert_type"

siginfo <- read.delim(file = "LDS-1613/Data/siginfo.txt")
colnames(siginfo)
# [1] "sig_id"                       "pert_id"                      "cmap_name"                   
# [4] "bead_batch"                   "nearest_dose"                 "pert_dose"                   
# [7] "pert_dose_unit"               "pert_idose"                   "pert_itime"                  
# [10] "pert_time"                    "pert_time_unit"               "cell_mfc_name"               
# [13] "pert_mfc_id"                  "nsample"                      "cc_q75"                      
# [16] "ss_ngene"                     "tas"                          "pct_self_rank_q25"           
# [19] "wt"                           "median_recall_rank_spearman"  "median_recall_rank_wtcs_50"  
# [22] "median_recall_score_spearman" "median_recall_score_wtcs_50"  "batch_effect_tstat"          
# [25] "batch_effect_tstat_pct"       "is_hiq"                       "qc_pass"                     
# [28] "pert_type"                    "cell_iname"                   "det_wells"                   
# [31] "det_plates"                   "distil_ids"                   "build_name"                  
# [34] "project_code"                 "is_exemplar_sig"              "is_ncs_sig"                  
# [37] "is_null_sig"          

instinfo <- read.delim(file = "LDS-1613/Data/instinfo.txt")
colnames(instinfo)
# [1] "sample_id"      "pert_id"        "cmap_name"      "bead_batch"     "nearest_dose"   "pert_dose"     
# [7] "pert_dose_unit" "pert_idose"     "pert_time"      "pert_itime"     "pert_time_unit" "cell_mfc_name" 
# [13] "pert_mfc_id"    "det_plate"      "det_well"       "rna_plate"      "rna_well"       "count_mean"    
# [19] "count_cv"       "qc_f_logp"      "qc_iqr"         "qc_slope"       "pert_type"      "cell_iname"    
# [25] "qc_pass"        "dyn_range"      "inv_level_10"   "build_name"     "failure_mode"   "project_code" 

# load gctx files
##------------------------------------------------------------------------------

# gctx
ds_path <- "LDS-1613/Data/level5_beta_all_n1102210x12328.gctx"
col_ids <- read_gctx_ids(gctx_path = ds_path, dim = "col")
row_ids <- read_gctx_ids(ds_path, dim = "row")

# load one drug at a time...
triptolide_meta <- siginfo[siginfo$cmap_name == "triptolide",]
triptolide_meta <- triptolide_meta[which(triptolide_meta$pert_time == "24"),]
triptolide_columns <- triptolide_meta$sig_id
col_meta <- read_gctx_meta(ds_path, dim = "col")
idx <- which(col_meta$id%in%triptolide_columns)
triptolide_ds <- parse_gctx(ds_path, cid = idx)
triptolide_mat <- as.data.frame(mat(triptolide_ds))

# found gene id conversion list on harmonizome...
# need to filter for L1000 measured genes...

# geneName mapping
##------------------------------------------------------------------------------
gene_list_terms <- read.table(file = "LDS-1613/Data/gene_list_terms.txt", header = T)[,c(1,3)]
#head(gene_list_terms)
# GeneSym GeneID
# 1    A1CF  29974
# 2     A2M      2
# 3  A4GALT  53947
# 4   A4GNT  51146
# 5    AAAS   8086
# 6    AACS  65985

# load previous release signature data to get genes to filter on...
phase_1_consensus <- 
  read.table(file = "LDS-1613/Data/matPH3_2_1_0.2_0.3_L1000_Batch2017_Regina_removed.txt", 
             header = T, row.names = 1)
L1000_measured <- colnames(phase_1_consensus)
L1000_measured_terms <- gene_list_terms[which(gene_list_terms$GeneSym %in% L1000_measured),]

measured_mat <- triptolide_mat[which(rownames(triptolide_mat) %in% L1000_measured_terms$GeneID),]
toRename <- as.data.frame(rownames(measured_mat))
colnames(toRename) <- c("GeneID")
matched <- merge(x = toRename, y = L1000_measured_terms, by.x = "GeneID", by.y = "GeneID", all.x = T)
rownames(measured_mat) <- matched$GeneSym

# overview heatmap
##------------------------------------------------------------------------------
colnames(triptolide_meta)
triptolide_anno_columns <- c("sig_id",                    
                               "bead_batch",                  
                               "pert_idose", 
                               "pert_itime",                  
                               "pert_time", "cell_iname")
triptolide_anno <- triptolide_meta[triptolide_anno_columns]
rownames(triptolide_anno) <- triptolide_anno$sig_id
triptolide_anno$sig_id <- NULL
pheatmap::pheatmap(measured_mat, fontsize = 4, annotation_col = triptolide_anno, show_rownames = F, show_colnames = F)

barplot(rowMeans(measured_mat))

################################################################################

# Filtering

# Calculating consensus signatures

#> The L1000 Level 4 datasets was filtered for only 24-hour treatment samples. 

# filtered prior to gctx load
# time_24h <- triptolide_meta[which(triptolide_meta$pert_time == "24"),]

#> Gene expression profiles were aggregated for samples using both the same small 
#> molecule and the same cell line (technical/biological replicates and/or use 
#> different doses of the small molecule). Aggregation was performed by first 
#> counting the number of samples that have a |z-score| > 1 for a particular gene. 

# aggregate cell line data
##------------------------------------------------------------------------------

# transpose matrix
t_measured <- t(measured_mat)
t_measured

# merge w/ pert_dose
tomerge <- as.data.frame(triptolide_meta$cell_iname, row.names = triptolide_meta$sig_id)
colnames(tomerge) <- c("cell_iname")

t_merge <- merge(t_measured, tomerge, by.x = "row.names", by.y = "row.names")
rownames(t_merge) <- t_merge$Row.names
t_merge$Row.names <- NULL

t_merge$cell_iname
length(colnames(t_merge))
colnames(t_merge)[which(colnames(t_merge) != "cell_iname")]
length(colnames(t_merge)[which(colnames(t_merge) != "cell_iname")])
# cell_agg <- t_merge %>% group_by("cell_iname") %>% mutate(across(c(colnames(t_merge)[which(colnames(t_merge) != "cell_iname")]), .fns = mean)) %>% distinct()
dim(t_merge)
t_merge$cell_iname

# cell_agg <- t_merge %>% group_by("cell_iname") %>% mutate(.fns = "mean") %>% distinct()

#all columns not cell_iname
library(Hmisc)
length(colnames(t_merge))
cols <- colnames(t_merge)[which(colnames(t_merge) %nin% c("cell_iname"))] 
length(cols)
cell_agg <- t_merge %>% group_by(cell_iname) %>% summarise(across(all_of(cols), mean), .groups = 'drop') %>% distinct() %>% as.data.frame()

dim(cell_agg)

rownames(cell_agg) <- cell_agg$cell_iname
cell_agg$cell_iname <- NULL

cell_agg
pheatmap::pheatmap(cell_agg, show_colnames = F)

#
##------------------------------------------------------------------------------
#> If this count was more than 20% of the total number of samples, then that 
#> particular gene was included in the aggregated expression profile.

#> Next, the above-aggregated gene expression profiles are collapsed at the 
#> small molecule level. Gene expression profiles that correspond to the same 
#> small molecule across all cell lines are aggregated to produce the TCS. A 
#> gene is included in the final TCS if it is up/down-regulated by a |z-score|
#>  > 1 in more than 30% of the cell lines that were treated by the same small
#>   molecule.

# generate list of genes to keep with |z-score| > 1 in 30% of cell lines

# how many total cell lines?
total_cell_lines <- length(rownames(cell_agg)) # 102

t_cell_agg <- as.data.frame(t(cell_agg))

pheatmap::pheatmap(t_cell_agg)


# keep_genes <- list()
# for (i in 1:length(rownames(t_cell_agg))){
#   # print(rownames(t_cell_agg)[i])
#   pass_filter_count <- list()
#   for (j in 1:length(colnames(t_cell_agg))){
#     # print(colnames(t_cell_agg)[j])
#     test_df <- t_cell_agg[i,j]
#     # head(test_df)
#     if (abs(test_df) > 1) {
#       # pass_filter_count <- paste(pass_filter_count, rownames(t_cell_agg)[i], sep="")
#       pass_filter_count <- append(pass_filter_count, rownames(t_cell_agg)[i])
#     }
#   }
#   if (length(pass_filter_count) > 0.3*total_cell_lines){
#     # keep_genes <- paste(keep_genes, rownames(t_cell_agg)[i], sep="")
#     keep_genes <- append(keep_genes, rownames(t_cell_agg)[i])
#   }
# }
# 
# length(keep_genes)

# Calculate the threshold for the number of columns
threshold <- ceiling(0.3 * ncol(t_cell_agg))

# Initialize an empty vector to store row indices that meet the threshold
selected_rows <- c()

# Loop through the rows
for (i in 1:nrow(t_cell_agg)) {
  # Count the number of columns with absolute values greater than 1
  count <- sum(abs(t_cell_agg[i, ]) > 1)
  
  # Check if the count meets the threshold
  if (count >= threshold) {
    selected_rows <- c(selected_rows, i)
  }
}

# Print the selected rows
print(selected_rows)

pheatmap::pheatmap(t_cell_agg)
pheatmap::pheatmap(t_cell_agg[selected_rows,])

barplot(rowMeans(t_cell_agg[selected_rows,]))

old_triptolide <- phase_1_consensus[c("triptolide"),]
old_triptolide_nonZero <- old_triptolide[, colSums(old_triptolide != 0) > 0]
dim(old_triptolide_nonZero)



rowMeans <- rowMeans(t_cell_agg)
threshold <- 1
rows_to_label <- which(rowMeans > threshold)

order_rows <- order(rowMeans, decreasing = T)
order_rows

length(rows_to_label)

pheatmap::pheatmap(t_cell_agg[order_rows,], cellwidth = 10,
                   cluster_rows = F,
                   # clustering_distance_rows = "euclidean",
                   # clustering_method = "complete",
                   # annotation_row = data.frame(Label = ifelse(1:nrow(t_cell_agg) %in% rows_to_label, "Above Threshold", "Below Threshold")),
                   show_rownames = FALSE  # Optional: Set to FALSE if you don't want the default row names
)

t_cell_agg <- t_cell_agg[which(rownames(t_cell_agg) %in% colnames(old_triptolide_nonZero)),]


mat_data <- as.matrix(t_cell_agg)
mat_data <- na.omit(mat_data)
row_means <- rowMeans(mat_data)
order_rows <- order(row_means, decreasing = T)
library(ComplexHeatmap)

dim(mat_data)
color_gradient <- colorRampPalette(c("blue", "white", "red"))(length(row_means))

pdf(file = "triptolide/triptolide_TCS_heatmap.pdf")
ComplexHeatmap::Heatmap(mat_data[order_rows,], width = unit(10, "cm"), show_row_names = F,
                        cluster_rows = F, right_annotation = ComplexHeatmap::rowAnnotation(TCS = anno_barplot(x = row_means[order_rows], fill = color_gradient, bar_width = 1)))
dev.off()

write.csv(mat_data, file = "triptolide/triptolide_l1000_24hour_processed.csv")


################################################################################
# Explore triptolide


# ds_path <- c("/data/ccbb/L1000/level5_beta_trt_cp_n720216x12328.gctx")
# 
# # first parse metadata to ID where relevant compounds are...
# col_meta <- read_gctx_meta(ds_path, dim="col")
# row_meta <- read_gctx_meta(ds_path, dim = "row")
# # figure out which signatures correspond to vorinostat by searching the 'pert_iname' column
# idx <- which(col_meta$pert_iname=="Zotiraciclib")
# 
# # read only those columns from the GCTX file by using the 'cid' parameter
# zotiraciclib_ds <- parse_gctx(ds_path, cid=idx)
