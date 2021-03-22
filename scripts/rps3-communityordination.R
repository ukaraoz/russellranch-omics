library(easypackages)
packages <- c("dplyr", "tibble", "magrittr", "tidyr", "labdsv", "ggplot2", "phyloseq", 
              "ecodist", "doParallel", "vegan", "DESeq2")
libraries(packages)
base = "/Users/ukaraoz/Work/russellranch"

source(file.path(base, "scripts/misc.R"))
outdir = file.path(base, "results")

metadata = read_metadata(file.path(base, "data/sample2metadata.txt"))
rownames(metadata) = metadata[,"sample"]

#tree.file = file.path(base, "msa_phylotree/RussellRanch_ALL.rps3_centroids.muscle.95gapsremoved.fasttree.tree")
tree.file = file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3_centroids.muscle.95gapsremoved.iqtree.tree")

rps3centroid.file = file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3centroid2longestscaffold.txt")
rps3centroid = read.table(rps3centroid.file, header = T, sep = "\t") %>% as_tibble

assembly_stats = readRDS(file.path(base, "data/assembly/assembly_stats.rds")) %>% 
  as.data.frame %>% rownames_to_column(var = "sample") %>% as_tibble

coverage.file = file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3meancoverage_99pid.txt")
coverage = read.table(coverage.file, check.names = F, header = T, sep = "\t") %>% tbl_df
colnames(coverage) = gsub("RussellRanch_ALL.rps3_centroidslongestscaffold.vs.|.sorted Mean", "", colnames(coverage))
coverage_norm = coverage
for(c in 2:ncol(coverage_norm)) {
  nreads = assembly_stats %>% filter(sample == colnames(coverage)[2]) %>% select(reads_prefiltered) %>% as.numeric
  coverage_norm[,c] = coverage[, c]/nreads*100000000
}
coverage_norm = coverage_norm %>%
  inner_join(rps3centroid, by = c("Contig" = "cluster.longestscaffold")) %>%
  select(c("Contig", "centroid.orf", everything())) %>%
  select(-"cluster", -"cluster.size", -"cluster.longestscaffoldlength")

# 1652 SGs starting

#SG coverage values were Hellinger standardized, and then SGs
#were removed that had a coefficient of variation (CV) of normalized coverage
#>3 or with <5 samples where raw coverage was ≥0.25

coverage_norm_matrix = coverage_norm %>% select(metadata[, "sample"]) %>% as.matrix
rownames(coverage_norm_matrix) = coverage_norm %>% select(`centroid.orf`) %>% pull
rawcov_threshold = 0.25
samplesover_threshold = 6
cv_threshold = 3
coverage_norm_matrix_thresholded = binarize(coverage_norm_matrix, rawcov_threshold)
centroid.orf.passing = rownames(coverage_norm_matrix_thresholded)[which(apply(coverage_norm_matrix_thresholded, 1, sum) > samplesover_threshold)]
# 1085 passing
coverage_norm_matrix = coverage_norm_matrix[centroid.orf.passing,]
coverage_norm_matrix_hellinger = decostand(coverage_norm_matrix, method = "hellinger")
coverage_norm_matrix_hellinger_cv = apply(coverage_norm_matrix_hellinger, 1, cv)
centroid.orf.passing = names(coverage_norm_matrix_hellinger_cv)[which(coverage_norm_matrix_hellinger_cv < cv_threshold)]
# 1034 passing
coverage_norm_matrix_hellinger_passing = coverage_norm_matrix_hellinger[centroid.orf.passing, ]

otutable = otu_table(coverage_norm_matrix_hellinger_passing, taxa_are_rows=TRUE)
tree = read_tree(tree.file)
phyloseqObject = phyloseq(otutable, tree, sample_data(metadata))

registerDoParallel(makeCluster(8))
all.udist.w = phyloseq::distance(phyloseqObject,"unifrac",weighted=TRUE)
wUnifracDist = as.matrix(all.udist.w)
sol = metaMDS(wUnifracDist, k = 2, try = 500, trymax = 500)
cat("Stress=", sol$stress, "\n") # Stress= 0.1043745
plotData = sol$points[, 1:2]

pdf(file = file.path(outdir, paste("nmdsonUnifracW.pdf", sep = "")), width = 12, height = 12)
par(mfrow=c(2,2), mar = c(2, 2, 2, 2))
colorbylist = c("color_by_depth", "color_by_treatment", "color_by_year", "color_by_plot")
for(i in 1:length(colorbylist)) {
  colorby = colorbylist[i]
  #pdf(file = file.path(outdir, paste("nmdsonUnifracW_", colorby, ".pdf", sep = "")), width = 8, height = 8)
  colors = as.matrix(sample_data(phyloseqObject)[rownames(as.matrix(wUnifracDist)), colorby])[,1]
  pch = as.matrix(sample_data(phyloseqObject)[rownames(as.matrix(wUnifracDist)), "pch_by_year"])[,1]
  plot(plotData, col = colors, pch = pch, xlab = "NMDS1", ylab = "NMDS2", cex = 1.8)
  temp = unique(metadata[, c(sub("color_by_", "", colorby), colorby)])
  temp1 = unique(metadata[, c("year", "pch_by_year")])

  legend("topright", cex = 1.4, 
          pch = 16,
         legend = temp[,1],
         col = temp[,2])
  legend("bottomright", cex = 1.8, 
          pch = temp1[, "pch_by_year"],
         legend = temp1[, "year"],
         col = "black")
  #
}
dev.off()

for(i in 1:length(colorbylist)) {
  mrpp <- with(sample_data(phyloseqObject), mrpp(wUnifracDist, get(sub("color_by_", "", colorbylist[i])), permutations = 10000, weight.type = 1))
  cat("mrpp-", sub("color_by_", "", colorbylist[i]), "A=", mrpp$A, " ", "p=", mrpp$Pvalue, "\n")

  # phyloseqObject_subset1 = subset_samples(phyloseqObject, depth %in% c("0-15", "15-30"))
  # samples_subset1 = sample_data(phyloseqObject_subset1) %>% pull(sample)
  # mrpp <- with(sample_data(phyloseqObject_subset1), 
  #            mrpp(wUnifracDist[samples_subset1,samples_subset1], get(sub("color_by_", "", colorbylist[i])), permutations = 10000, weight.type = 1))
}
# mrpp- depth A= 0.2059556   p= 9.999e-05 
# mrpp- treatment A= 0.0007550516   p= 0.3784622 
# mrpp- year A= 0.03353623   p= 0.00129987 
# mrpp- plot A= 0.02383107   p= 0.1643836 

#nmds <- nmds(wUnifracDist.DNA, mindim=2, maxdim=2)
#nmin <- nmds.min(nmds)
# Minimum stress for given dimensionality:  0.1545646 
# r^2 for minimum stress configuration:  0.9203968 
#plotData = cbind(nmin[[1]], nmin[[2]])


