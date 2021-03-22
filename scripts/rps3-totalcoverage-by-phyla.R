library(easypackages)
packages <- c("dplyr", "tibble", "magrittr", "tidyr", "labdsv", "ggplot2", "phyloseq", "stringr",
              "ecodist", "doParallel", "vegan", "DESeq2")
libraries(packages)
base = "/Users/ukaraoz/Work/russellranch"

source(file.path(base, "scripts/misc.R"))
metadata = read_metadata(file.path(base, "data/sample2metadata.txt"))

depths = levels(metadata$depth)
treatment_depths = levels(metadata$treatment_depth)

dereplicated_genomes_comp50cont25 = read.table(file.path(base, "data/binning/dereplicated_genomes_comp50cont25.txt")) %>% as_tibble
colnames(dereplicated_genomes_comp50cont25) = "bin"
dereplicated_genomes_comp50cont25 = dereplicated_genomes_comp50cont25 %>%
  mutate(`bin` = str_replace(`bin`, ".fa", ".fa.faa"))

#rp16_gene2orf2bin = readRDS(file.path("/Users/ukaraoz/Work/russellranch/results/annotation/phylogeny_rp16/", "rp16_gene2orf2bin.rds"))

rps3centroid = read.table(file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3centroid2longestscaffold.txt"), header = T, sep = "\t") %>% as_tibble
rps3_centroids2scaffold2length = read.table(file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3_centroids2scaffold2length.txt"), 
                                            header = T, sep = "\t") %>% as_tibble
rps3_centroids2othermembers = rps3_centroids2scaffold2length %>% 
  filter(is_centroid == TRUE) %>%
  inner_join(rps3_centroids2scaffold2length, by = c("cluster" = "cluster")) %>%
  #filter(is_centroid.x == TRUE & is_centroid.y == FALSE) %>%
  select(cluster, orf.x, orf.y) %>%
  dplyr::rename(centroid.orf = `orf.x`, linked.orf = `orf.y`)

rps3_2_taxa = read.table(file.path(base, "data/annotation/phylogenetic/rps3/rps3_2_taxa.txt"), header = T, sep = "\t") %>% as_tibble
assembly_stats = readRDS(file.path(base, "data/assembly/assembly_stats.rds")) %>% 
  as.data.frame %>% rownames_to_column(var = "sample") %>% as_tibble

coverage = read.table(file.path(base, "data/annotation/phylogenetic/rps3/RussellRanch_ALL.rps3meancoverage_99pid.txt"), check.names = F, header = T, sep = "\t") %>% as_tibble
colnames(coverage) = gsub("RussellRanch_ALL.rps3_centroidslongestscaffold.vs.|.sorted Mean", "", colnames(coverage))
coverage_norm = coverage
for(c in 2:ncol(coverage_norm)) {
  nreads = assembly_stats %>% filter(sample == colnames(coverage)[2]) %>% select(reads_prefiltered) %>% as.numeric
  coverage_norm[,c] = coverage[, c]/nreads*100000000
}

coverage_norm_long = coverage_norm %>% 
  gather("sample", "coverage", -`Contig`, factor_key=FALSE) %>%
  inner_join(metadata, by = c("sample" = "sample")) %>%
  select(c("Contig", "sample", "coverage", "depth", "treatment_depth"))

# total coverage
coverage_norm_totalcov = coverage_norm_long %>% 
  group_by(`Contig`) %>%
  summarise(`total coverage` = sum(coverage))

# by depth
coverage_norm_long_bydepth = coverage_norm_long %>% 
  select(-treatment_depth, -sample) %>%
  group_split(`depth`)
coverage_norm_bydepth = do.call("rbind", coverage_norm_long_bydepth) %>% tbl_df %>%
  group_by(`Contig`, `depth`) %>%
  summarise(`coverage` = sum(coverage)) %>%
  ungroup() %>% as.data.frame %>%
  tidyr::pivot_wider(names_from = depth, values_from = coverage)

# by depth and treatment
coverage_norm_long_bytreatment_depth = coverage_norm_long %>% 
  select(-depth, -sample) %>%
  group_split(`treatment_depth`)
coverage_norm_bytreatment_depth = do.call("rbind", coverage_norm_long_bytreatment_depth) %>% tbl_df %>%
  group_by(`Contig`, `treatment_depth`) %>%
  summarise(`coverage` = sum(coverage)) %>%
  ungroup() %>% as.data.frame %>%
  tidyr::pivot_wider(names_from = treatment_depth, values_from = coverage)

percent <- function(x) (x / sum(x) * 100)
coverage_norm_byfactors = coverage_norm_totalcov %>%
  inner_join(coverage_norm_bydepth, by = c("Contig" = "Contig")) %>%
  inner_join(coverage_norm_bytreatment_depth, by = c("Contig" = "Contig")) # %>% mutate_at(vars(-Contig), percent)

#saveRDS(coverage_norm_percentcovbyfactors, 
#         file.path("/Users/ukaraoz/Work/russellranch/results/annot/rps3/rps3.coverage_norm_percentcovbyfactors.rds"))

#coverage_norm_percentcovbyfactors = 
#  readRDS("/Users/ukaraoz/Work/russellranch/results/annot/rps3/rps3.coverage_norm_percentcovbyfactors.rds")

coverage_norm_byfactors = coverage_norm_byfactors %>%
  inner_join(rps3centroid, by = c("Contig" = "cluster.longestscaffold")) %>%
  #select(c("Contig", "total coverage", "centroid.orf")) %>%
  inner_join(rps3_2_taxa, by = c("centroid.orf" = "rps3.orf"))

coverage_norm_byfactors_byphyla = coverage_norm_byfactors %>%
  group_by(`taxa`) %>%
  summarise(`total coverage` = sum(`total coverage`),
            `0-15` = sum(`0-15`),
            `15-30` = sum(`15-30`),
            `30-60` = sum(`30-60`),
            `60-100` = sum(`60-100`))

threshold_75percent_total = quantile(coverage_norm_byfactors %>% pull(`total coverage`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["75%"]
threshold_90percent_total = quantile(coverage_norm_byfactors %>% pull(`total coverage`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["90%"]

threshold_75percent_0_15 = quantile(coverage_norm_byfactors %>% pull(`0-15`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["75%"]
threshold_90percent_0_15 = quantile(coverage_norm_byfactors %>% pull(`0-15`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["90%"]

threshold_75percent_15_30 = quantile(coverage_norm_byfactors %>% pull(`15-30`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["75%"]
threshold_90percent_15_30 = quantile(coverage_norm_byfactors %>% pull(`15-30`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["90%"]

threshold_75percent_30_60 = quantile(coverage_norm_byfactors %>% pull(`30-60`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["75%"]
threshold_90percent_30_60 = quantile(coverage_norm_byfactors %>% pull(`30-60`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["90%"]

threshold_75percent_60_100 = quantile(coverage_norm_byfactors %>% pull(`60-100`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["75%"]
threshold_90percent_60_100 = quantile(coverage_norm_byfactors %>% pull(`60-100`), probs = seq(0, 1, 0.05), na.rm = FALSE, names = TRUE)["90%"]

coverage_norm_byfactors = coverage_norm_byfactors %>%
  mutate(greaterthan_threshold_75percent_total = ifelse(`total coverage` >= threshold_75percent_total, "red", "gray")) %>%
  mutate(greaterthan_threshold_90percent_total = ifelse(`total coverage` >= threshold_90percent_total, "red", "gray")) %>%
  mutate(greaterthan_threshold_75percent_0_15 = ifelse(`0-15` >= threshold_75percent_0_15, "red", "gray")) %>%
  mutate(greaterthan_threshold_90percent_0_15 = ifelse(`0-15` >= threshold_90percent_0_15, "red", "gray")) %>%
  mutate(greaterthan_threshold_75percent_15_30 = ifelse(`15-30` >= threshold_75percent_15_30, "red", "gray")) %>%
  mutate(greaterthan_threshold_90percent_15_30 = ifelse(`15-30` >= threshold_90percent_15_30, "red", "gray")) %>%
  mutate(greaterthan_threshold_75percent_30_60 = ifelse(`30-60` >= threshold_75percent_30_60, "red", "gray")) %>%
  mutate(greaterthan_threshold_90percent_30_60 = ifelse(`30-60` >= threshold_90percent_30_60, "red", "gray")) %>%
  mutate(greaterthan_threshold_75percent_60_100 = ifelse(`60-100` >= threshold_75percent_60_100, "red", "gray")) %>%
  mutate(greaterthan_threshold_90percent_60_100 = ifelse(`60-100` >= threshold_90percent_60_100, "red", "gray")) %>%
  inner_join(coverage_norm_byfactors_byphyla, by = c("taxa" = "taxa"))

namedcolors = c("gray", "red")
names(namedcolors) = c("gray", "red")
ggplot(coverage_norm_byfactors) + 
    geom_bar(aes(x=`percent total coverage`, y = reorder(taxa, `taxa total coverage`), 
                                   fill = `greaterthan_threshold_90percent`, 
                                   group = -desc(`percent total coverage`)),
             position=position_stack(reverse = FALSE), stat="identity", color="black", size=0.05) +
    scale_fill_manual(values = namedcolors) +
    xlab("percent total coverage") +
    ylab("taxa") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          axis.title.x = element_text(size = 22),
          axis.title.y = element_text(size = 22),
          legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(6, 6, 6, 6))


# = dereplicated_genomes_comp50cont10 %>%
# inner_join(rp16_gene2orf2bin, by = c("bin" = "bin")) %>%
#overage_norm_totalcov %>%
# inner_join(rp16_gene2orf2bin,) %>% 


# inner_join(dereplicated_genomes_comp50cont10, by = c("tree.bin" = "bin"))  %>% # 171
# inner_join(rp16_gene2orf2bin, by = c("tree.bin" = "bin")) %>%
# filter(gene == "rpsC") %>% # 152  # there are a few 
# inner_join(rps3_centroids2othermembers, by = c("gene_name" = "linked.orf")) %>%



