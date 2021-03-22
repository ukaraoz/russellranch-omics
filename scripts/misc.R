cv <- function(x, na.rm = FALSE)  {
  sd(x, na.rm = na.rm)/mean(x, na.rm = na.rm)
}

binarize=function(x,threshold=NA) {
  matd=x
  if(is.na(threshold))
    {
    threshold=min(x)+(max(x)-min(x))/2
    print(paste("Threshold: ",threshold))
    }

  matd[matd<=threshold]=0
  matd[matd>threshold]=1
  matd
}

read_metadata <- function(file) {
  sample2metadata = read.table(file, header = T, check.names = F, stringsAsFactors = F)
  sample2metadata = sample2metadata[order(as.numeric(sub("S", "", sample2metadata[, "sampleid"]))),]
  metadata = data.frame(sample = as.character(sample2metadata[, "sample"]),
                        sampleid = as.character(sample2metadata[, "sampleid"]),
                        replicate = factor(sample2metadata[, "replicate"],
                                           levels = c("A", "B", "C"),
                                           ordered = F),
                        plot = factor(sample2metadata[, "plot"],
                                      levels = c("1-1", "1-2", "1-4", "2-3", "2-4", 
                                                 "4-5", "5-5", "6-3", "6-4", "6-5", 
                                                 "6-7", "6-8", "6-9", "7-8", "8-8", 
                                                 "8-9"),
                                      ordered = F),
                        treatment = factor(as.character(sample2metadata[, "treatment"]),
                                                 levels = c("CMT", "LMT", "OMT"),
                                                 ordered = F),
                        depth = factor(as.character(sample2metadata[, "depth"]),
                                                 levels = c("0-15", "15-30", "30-60", "60-100"),
                                                 ordered = F),
                        year = factor(as.character(sample2metadata[, "year"]),
                                                       levels = c("Even", "Odd"),
                                                       ordered = F),
                        treatment_depth = factor(as.character(sample2metadata[, "treatment_depth"]),
                                                       ordered = F),
                        treatment_depth_year = factor(as.character(sample2metadata[, "treatment_depth_year"]),
                                                       ordered = F),
                        color_by_depth = as.character(sample2metadata[, "color_by_depth"]),
                        color_by_treatment = as.character(sample2metadata[, "color_by_treatment"]),
                        color_by_year = as.character(sample2metadata[, "color_by_year"]),
                        color_by_plot = as.character(sample2metadata[, "color_by_plot"]),
                        pch_by_treatment = as.character(sample2metadata[, "pch_by_treatment"]),
                        pch_by_year = as.numeric(sample2metadata[, "pch_by_year"]),
                        stringsAsFactors = F)
  return(metadata)
}

read_coverm_readcounts <- function(file) {
  #file = "/Users/ukaraoz/Work/russellranch/results/annot/rps3/meancoverage_99pid.readcounts.txt"
  readcounts = read.table(file, sep = "\t", header = T)
  colnames(readcounts) = gsub("RussellRanch_ALL.rps3_centroidslongestscaffold.vs.|.sorted.Read.Count", "", colnames(readcounts))
  readcounts = readcounts %>% tbl_df %>%
    dplyr::rename(cluster.longestscaffold = Contig) %>%
    dplyr::inner_join(rps3centroid, by = c("cluster.longestscaffold" = "cluster.longestscaffold")) %>%
    dplyr::select(-cluster, -cluster.size, -cluster.longestscaffoldlength) %>%
    dplyr::select(c("cluster.longestscaffold", "centroid.orf", everything()))
  return(readcounts)
}

p2signiflevels <- function(p) {
  map_signif_level <- c("****"=0.0001, "***"=0.001, "**"=0.01,  "*"=0.05, "NS"=1)
  if(is.na(p)) {
    signif_level = 1
  } else {
    signif_level <- names(which.min(map_signif_level[which(map_signif_level >= p)]))
  }
  signif_level
}

estimate2sign <- function(a) {
  if(is.na(a)) { 
    sign = "0" 
  } else {
    if(a > 0) { sign = "+" }
    if(a < 0) { sign = "-" }
    if(a == 0) { sign = "0" }
  }
  sign
}
