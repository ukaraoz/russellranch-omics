library(easypackages)
packages <- c("dplyr", "tibble", "magrittr", "tidyr", "labdsv", "ggplot2", "phyloseq", 
              "ecodist", "doParallel", "vegan", "DESeq2")
libraries(packages)
base = "/Users/ukaraoz/Work/russellranch"

source(file.path(base, "scripts/misc.R"))
outdir = file.path(base, "results")
verbose = 0

rps3centroid.file = file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3centroid2longestscaffold.txt")
rps3centroid = read.table(rps3centroid.file, header = T, sep = "\t") %>% as_tibble

metadata = read_metadata(file.path(base, "data/sample2metadata.txt"))


rps3centroid_readcounts = readRDS(file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3centroid_instrain_readcounts.rds"))
# rps3centroid_coverm_readcounts = read_coverm_readcounts("/Users/ukaraoz/Work/russellranch/data/annotation/rps3/meancoverage_99pid.readcounts.txt")
rps3centroid_readcounts_matrix = rps3centroid_readcounts %>% 
  select(metadata[, "sample"]) %>% as.matrix
rownames(rps3centroid_readcounts_matrix) = rps3centroid_readcounts %>% select(`centroid.orf`) %>% pull

### Build DEseq Object - Depth
dds <- DESeqDataSetFromMatrix(countData = rps3centroid_readcounts_matrix,
                              colData = metadata,
                              design = ~ replicate + treatment + year + depth)
# print(dds)

# run DESeq with LRT
dds_lrt_out <- DESeq(dds, test = "LRT", 
                     reduced = ~ replicate + treatment + year)
# Output Size Factors
s_factors_lrt <- sizeFactors(dds_lrt_out)

write.table(s_factors_lrt, file.path(outdir, "rps3-enrichment-depth_size_factors_lrt.txt"), col.names =F, sep = "\t",quote = FALSE)
# Plot Dispersion Estimates
pdf(file.path(outdir, "rps3-enrichment-depth_dispersionestimates.pdf"))
plotDispEsts(dds_lrt_out)
dev.off()

### Output Primary Results Table
# Show Results Variables
resultsNames(dds_lrt_out)

# Display Differential Abundance Summary
summary(results(dds_lrt_out, alpha = 0.05))

# Write Filtered Results (FDR < 0.05) To Disk
#depth_LRT_results_p <- subset(results(dds_lrt_out, alpha = 0.05, format = "DataFrame"), padj < 0.05)
depth_LRT_results <- results(dds_lrt_out, alpha = 0.05, format = "DataFrame")
write.table(as.data.frame(depth_LRT_results), file.path(outdir, "rps3-enrichment-depth_LRT_results_p.txt"), sep = "\t", quote = FALSE)

depth_LRT_results_p <- subset(results(dds_lrt_out, alpha = 0.05, format = "DataFrame"), padj < 0.05)

### Run Linear Model on All Significant Results to Check Slope Across Depth

# Convert Results to DataFrame and add Colums for Slopes 
results_table <- as.data.frame(depth_LRT_results)
results_table <- data.frame(results_table, lm_slope = NA, slope_p = NA)

# Loop Through all Significant Results 
for (i in 1:nrow(depth_LRT_results_p)) {
  cat("Testing ", i, " out of ", nrow(depth_LRT_results_p), " total", "\n")
  # Get Rowname of Row i
  tmp_OTU <- rownames(depth_LRT_results_p)[i]
  
  # Get Counts vs Depth Data for SG, Log Scale Counts, and Create Numeric Dummy Depth Variable
  tmp_counts <- plotCounts(dds, gene=tmp_OTU, intgroup=c("depth"),
                           returnData=TRUE)
  tmp_counts <- data.frame(tmp_counts,
                           log_count = log(tmp_counts$count), 
                           dummy_depth = as.numeric(tmp_counts$depth))
  
  # Run Linear Model log_count = m*dummy_depth + b
  lm_d_log_dum <- lm(log_count ~ dummy_depth, data = tmp_counts)
  tmp_summary <- summary(lm_d_log_dum)
  
  # Add Model Coefficients to DataFrame
  results_table[tmp_OTU, "lm_slope"] <- tmp_summary$coefficients["dummy_depth", "Estimate"]
  results_table[tmp_OTU, "slope_p"] <- tmp_summary$coefficients["dummy_depth", "Pr(>|t|)"] 
}
# Add SG Centroid Names to Results table
results_table = data.frame(centroid.orf = rownames(results_table),
                            results_table)
results_table_backup = results_table
# FDR Correct Slope p-values and Subset for FDR <= 0.05
results_table = results_table %>% as_tibble %>%
  mutate(slope_fdr = p.adjust(results_table$slope_p, n = length(which(results_table$padj<=0.05)), method = "fdr"))
results_table = results_table %>%
  mutate(depth_call_category = case_when(((results_table$slope_fdr > 0.05) | is.na(results_table$slope_fdr)) ~ "NS",
                                          (results_table$slope_fdr <= 0.05) & (results_table$lm_slope > 0) ~ "Increase",
                                          (results_table$slope_fdr <= 0.05) & (results_table$lm_slope < 0) ~ "Decrease"))
write.table(results_table, file.path(outdir, "rps3-enrichment-depth_LRT_results.xls"), row.names = F, col.names = T, sep = "\t", quote = FALSE)
saveRDS(results_table, file.path(outdir, "rps3-enrichment-depth_LRT_results.rds"))


####################
# 2. Differential Abundance by treatment
####################

#### Load Data -- Specific Combined Variables for Compairing Treatment at Each Depth
#colData  <- read.delim("~/Dropbox/Banfield_Lab_Files/Projects/Angelo2014/Studies/17_03_14_DEseq_of_m2_Mapping/Sample_Metadata_ALL_Grouped.txt",
#                       header = TRUE, sep = "\t") 

### Build DEseq Object - Treatment 
#dds <- DESeqDataSetFromMatrix(countData = countData,
#                              colData = colData,
#                              design = ~ Replicate + Time_Point + Factor)
#dds

depths = levels(metadata[, "depth"])
results = matrix(nrow = nrow(rps3centroid_readcounts_matrix), ncol = 0)
for(i in 1:length(depths)) {
  cat("Testing ", depths[i], "\n")
  depth.samples = subset(metadata, depth == depths[i])[, "sample"]
  dds <- DESeqDataSetFromMatrix(countData = rps3centroid_readcounts_matrix[, depth.samples],
                                colData = subset(metadata, sample %in% depth.samples),
                                design = ~ replicate + year + treatment)
  # print(dds)
  ### Run DEseq Stats with Standard Wald Test with Local Fitting
  dds_out <- DESeq(dds, fitType = "local")
  
  ### Output Size Factors and Plot Dispersion Estimates
  
  # Output Size Factors
  s_factors <- sizeFactors(dds_out)
  write.table(s_factors_lrt, 
              file.path(outdir, paste("rps3-enrichment-", "depth_", depths[i], "treatment_size_factors_wald.txt", sep = "")), 
              col.names =F, sep = "\t",quote = FALSE)
  
  # Plot Dispersion Estimates
  #pdf(file.path(outdir, "treatment_dispersionestimates.pdf")
  #plotDispEsts(dds_out)
  #dev.off()
  
  ### Output Results Tables
  
  # Show Results Variables
  resultsNames(dds_out)

  results_LMTvsCMT = as.data.frame(results(dds_out, contrast = c("treatment", "LMT", "CMT"), alpha = 0.05))
  colnames(results_LMTvsCMT) = paste(colnames(results_LMTvsCMT), ".", depths[i], ".", "LMTvsCMT", sep = "")
  results_OMTvsCMT = as.data.frame(results(dds_out, contrast = c("treatment", "OMT", "CMT"), alpha = 0.05))
  colnames(results_OMTvsCMT) = paste(colnames(results_OMTvsCMT), ".", depths[i], ".", "OMTvsCMT", sep = "")
  results = cbind(results,
                  results_LMTvsCMT, 
                  results_OMTvsCMT[rownames(results_LMTvsCMT), ])
}

results_table = data.frame(centroid.orf = rownames(results), check.names = F,
                           subset(results, select = grep("log2FoldChange|pvalue|padj", colnames(results))),
                           `log2FoldChange.0-15.LMTvsCMT.sign` = ifelse(is.na(results$`padj.0-15.LMTvsCMT`) | (results$`padj.0-15.LMTvsCMT` > 0.05), "NS", "S"),
                           `log2FoldChange.15-30.LMTvsCMT.sign` = ifelse(is.na(results$`padj.15-30.LMTvsCMT`) | (results$`padj.15-30.LMTvsCMT` > 0.05), "NS", "S"),
                           `log2FoldChange.30-60.LMTvsCMT.sign` = ifelse(is.na(results$`padj.30-60.LMTvsCMT`) | (results$`padj.30-60.LMTvsCMT` > 0.05), "NS", "S"),
                           `log2FoldChange.60-100.LMTvsCMT.sign` = ifelse(is.na(results$`padj.60-100.LMTvsCMT`) | (results$`padj.60-100.LMTvsCMT` > 0.05), "NS", "S"),
                           `log2FoldChange.0-15.OMTvsCMT.sign` = ifelse(is.na(results$`padj.0-15.OMTvsCMT`) | (results$`padj.0-15.OMTvsCMT` > 0.05), "NS", "S"),
                           `log2FoldChange.15-30.OMTvsCMT.sign` = ifelse(is.na(results$`padj.15-30.OMTvsCMT`) | (results$`padj.15-30.OMTvsCMT` > 0.05), "NS", "S"),
                           `log2FoldChange.30-60.OMTvsCMT.sign` = ifelse(is.na(results$`padj.30-60.OMTvsCMT`) | (results$`padj.30-60.OMTvsCMT` > 0.05), "NS", "S"),
                           `log2FoldChange.60-100.OMTvsCMT.sign` = ifelse(is.na(results$`padj.60-100.OMTvsCMT`) | (results$`padj.60-100.OMTvsCMT` > 0.05), "NS", "S")
                           ) %>% as_tibble
results_table = results_table %>%
  mutate(`log2FoldChange.0-15.LMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.0-15.LMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.0-15.LMTvsCMT.sign` == "S") & (`log2FoldChange.0-15.LMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.0-15.LMTvsCMT.sign` == "S") & (`log2FoldChange.0-15.LMTvsCMT` < 0) ~ "Decrease")) %>%
  mutate(`log2FoldChange.15-30.LMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.15-30.LMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.15-30.LMTvsCMT.sign` == "S") & (`log2FoldChange.15-30.LMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.15-30.LMTvsCMT.sign` == "S") & (`log2FoldChange.15-30.LMTvsCMT` < 0) ~ "Decrease")) %>%
  mutate(`log2FoldChange.30-60.LMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.30-60.LMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.30-60.LMTvsCMT.sign` == "S") & (`log2FoldChange.30-60.LMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.30-60.LMTvsCMT.sign` == "S") & (`log2FoldChange.30-60.LMTvsCMT` < 0) ~ "Decrease")) %>%
  mutate(`log2FoldChange.60-100.LMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.60-100.LMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.60-100.LMTvsCMT.sign` == "S") & (`log2FoldChange.60-100.LMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.60-100.LMTvsCMT.sign` == "S") & (`log2FoldChange.60-100.LMTvsCMT` < 0) ~ "Decrease")) %>%
  mutate(`log2FoldChange.0-15.OMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.0-15.OMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.0-15.OMTvsCMT.sign` == "S") & (`log2FoldChange.0-15.OMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.0-15.OMTvsCMT.sign` == "S") & (`log2FoldChange.0-15.OMTvsCMT` < 0) ~ "Decrease")) %>%
  mutate(`log2FoldChange.15-30.OMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.15-30.OMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.15-30.OMTvsCMT.sign` == "S") & (`log2FoldChange.15-30.OMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.15-30.OMTvsCMT.sign` == "S") & (`log2FoldChange.15-30.OMTvsCMT` < 0) ~ "Decrease")) %>%
  mutate(`log2FoldChange.30-60.OMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.30-60.OMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.30-60.OMTvsCMT.sign` == "S") & (`log2FoldChange.30-60.OMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.30-60.OMTvsCMT.sign` == "S") & (`log2FoldChange.30-60.OMTvsCMT` < 0) ~ "Decrease")) %>%
  mutate(`log2FoldChange.60-100.OMTvsCMT.sign.category` = 
     case_when(`log2FoldChange.60-100.OMTvsCMT.sign` == "NS" ~ "NS",
              (`log2FoldChange.60-100.OMTvsCMT.sign` == "S") & (`log2FoldChange.60-100.OMTvsCMT` > 0) ~ "Increase",
              (`log2FoldChange.60-100.OMTvsCMT.sign` == "S") & (`log2FoldChange.60-100.OMTvsCMT` < 0) ~ "Decrease")) %>%
  as.data.frame
write.table(results_table, file.path(outdir, "rps3-enrichment-treatment_wald_results.xls"), row.names = F, col.names = T, sep = "\t", quote = FALSE)
saveRDS(results_table, file.path(outdir, "rps3-enrichment-treatment_wald_results.rds"))
