library(easypackages)
packages <- c("dplyr", "tibble", "magrittr", "tidyr", "labdsv", "ggplot2", "phyloseq", 
              "ecodist", "doParallel", "vegan", "DESeq2", "grid", "gridExtra", "ggpubr")
libraries(packages)

base = "/Users/ukaraoz/Work/russellranch"
source(file.path(base, "scripts/misc.R"))
outdir = file.path(base, "results")

#outdir = "/Users/ukaraoz/Work/russellranch/results/annot/rps3/diffabundance"
sample2metadata = read_metadata(file.path(base, "data/sample2metadata.txt")) %>% as_tibble
depths = sample2metadata %>% dplyr::select("depth") %>% distinct %>% pull %>% as.character
treatments = sample2metadata %>% dplyr::select("treatment") %>% distinct %>% pull %>% as.character

# per depth: CMT=5, LMT=5, OMT=6
# sample2metadata %>% dplyr::filter(depth == depths[j]) %>% dplyr::select(treatment)

# scaffold2length = readRDS(file.path(base, "results/annot", "RussellRanch.assembly.slf.fa.length.rds"))
# scaffold2trait = readRDS(file.path(base, "results/annot/microtrait", "RussellRanch_all.microtrait.scaffold2trait.rds")) %>%
#   dplyr::inner_join(scaffold2length, by = c("scaffold" = "scaffold"))
# scaffoldover2K2trait = scaffold2trait %>%
#   dplyr::filter(length >= 2000)
# saveRDS(scaffoldover2K2trait, file.path(base, "results/annot/microtrait", "RussellRanch_all.microtrait.scaffoldover2K2trait.rds"))
scaffold2trait = readRDS(file.path(base, "data/annotation/functional/microtrait", "RussellRanch_all.microtrait.scaffoldover2K2trait.rds"))

scaffold2coverage = readRDS(file.path(base, "data/binning", "RussellRanch_all.scaffold2coverage.rds")) %>%
  dplyr::inner_join(sample2metadata, by = c("sample" = "sample")) %>%
  dplyr::select(c("contig", "sample", "coverage", "treatment", "depth", "year", "replicate", "plot")) %>%
  # compute relative abundances
  group_by(sample,depth) %>%
  mutate(coverage_perc= prop.table(coverage) * 100)

alltraits = scaffold2trait %>% dplyr::select(`trait-display-short`) %>% distinct %>% pull %>% as.character %>% sort

verbose = 0
tbl_colnames = c("trait", "depth", "comparison", "p", "estimate", "n_scaffoldswtrait", "n_scaffoldswtraitfromdepth")
treatment_sign = as_tibble(setNames(list(character(), character(), character(), numeric(), numeric(), integer(), integer()), tbl_colnames))
for(i in 1:length(alltraits)) {
  message(i, "\t", alltraits[i], "\n")
  scaffoldswtrait = scaffold2trait %>% 
    dplyr::filter(`trait-display-short` == alltraits[i]) %>% 
    dplyr::select(c("scaffold", "n", "trait-type", "length")) %>% 
    dplyr::mutate(`assembly-depth` = stringr::str_replace(scaffold, "(.*)_(.*)_(.*)_(.*)_(.*)_(.*)_(.*)$", "\\4"))

  if(verbose) { message(scaffoldswtrait %>% pull(scaffold) %>% unique %>% length, " scaffolds with trait ", alltraits[i]) }
  scaffoldswtrait = scaffoldswtrait %>% 
    dplyr::filter(length >= 2000)
  n_scaffoldswtrait = scaffoldswtrait %>% dplyr::select(scaffold) %>% n_distinct
  if(verbose) { message("\t", n_scaffoldswtrait, " scaffolds over 2K") }

  for(j in 1:length(depths)) {
    scaffoldswtraitfromdepth = scaffoldswtrait %>%
      dplyr::filter(`assembly-depth` == depths[j])
    n_scaffoldswtraitfromdepth = scaffoldswtraitfromdepth %>% dplyr::select(scaffold) %>% n_distinct
    if(verbose) { message("\t\t", scaffoldswtraitfromdepth %>% pull(scaffold) %>% unique %>% length, " scaffolds assembled from depth ", depths[j]) }
    
    #c = setdiff(b,a)[1:10]
    #scaffold2coverage %>% dplyr::filter(contig == c[1])
  
    if(n_scaffoldswtraitfromdepth != 0) {
      temp = scaffoldswtraitfromdepth %>%
        dplyr::inner_join(scaffold2coverage, by = c("scaffold" = "contig")) %>% 
        dplyr::filter(depth == depths[j]) %>%
        dplyr::select(c("coverage", "coverage_perc", "treatment", "year", "n", "trait-type")) %>%
        dplyr::mutate(coverage_weighted = `coverage`*`n`,
                      coverage_perc_weighted = `coverage_perc`*`n`)
      temp_CMT = temp %>% dplyr::filter(treatment == "CMT") %>% dplyr::select(coverage_weighted) %>% pull
      temp_LMT = temp %>% dplyr::filter(treatment == "LMT") %>% dplyr::select(coverage_weighted) %>% pull
      temp_OMT = temp %>% dplyr::filter(treatment == "OMT") %>% dplyr::select(coverage_weighted) %>% pull
      temp_CMT_perc = temp %>% dplyr::filter(treatment == "CMT") %>% dplyr::select(coverage_perc_weighted) %>% pull
      temp_LMT_perc = temp %>% dplyr::filter(treatment == "LMT") %>% dplyr::select(coverage_perc_weighted) %>% pull
      temp_OMT_perc = temp %>% dplyr::filter(treatment == "OMT") %>% dplyr::select(coverage_perc_weighted) %>% pull
      
      temp_CMT_odd = temp %>% dplyr::filter(treatment == "CMT" & year == "Odd") %>% dplyr::select(coverage_weighted) %>% pull
      temp_LMT_odd = temp %>% dplyr::filter(treatment == "LMT" & year == "Odd") %>% dplyr::select(coverage_weighted) %>% pull
      temp_OMT_odd = temp %>% dplyr::filter(treatment == "OMT" & year == "Odd") %>% dplyr::select(coverage_weighted) %>% pull
      temp_CMT_oddperc = temp %>% dplyr::filter(treatment == "CMT" & year == "Odd") %>% dplyr::select(coverage_perc_weighted) %>% pull
      temp_LMT_oddperc = temp %>% dplyr::filter(treatment == "LMT" & year == "Odd") %>% dplyr::select(coverage_perc_weighted) %>% pull
      temp_OMT_oddperc = temp %>% dplyr::filter(treatment == "OMT" & year == "Odd") %>% dplyr::select(coverage_perc_weighted) %>% pull

      temp_CMT_even = temp %>% dplyr::filter(treatment == "CMT" & year == "Even") %>% dplyr::select(coverage_weighted) %>% pull
      temp_LMT_even = temp %>% dplyr::filter(treatment == "LMT" & year == "Even") %>% dplyr::select(coverage_weighted) %>% pull
      temp_OMT_even = temp %>% dplyr::filter(treatment == "OMT" & year == "Even") %>% dplyr::select(coverage_weighted) %>% pull
      temp_CMT_evenperc = temp %>% dplyr::filter(treatment == "CMT" & year == "Even") %>% dplyr::select(coverage_perc_weighted) %>% pull
      temp_LMT_evenperc = temp %>% dplyr::filter(treatment == "LMT" & year == "Even") %>% dplyr::select(coverage_perc_weighted) %>% pull
      temp_OMT_evenperc = temp %>% dplyr::filter(treatment == "OMT" & year == "Even") %>% dplyr::select(coverage_perc_weighted) %>% pull

      # CMT vs LMT
      wilcox.test.temp = wilcox.test(temp_CMT, temp_LMT, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMT.vs.CMT", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)

      wilcox.test.temp = wilcox.test(temp_CMT_perc, temp_LMT_perc, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTperc.vs.CMTperc", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      ## by year
      wilcox.test.temp = wilcox.test(temp_CMT_odd, temp_LMT_odd, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTodd.vs.CMTodd", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      wilcox.test.temp = wilcox.test(temp_CMT_oddperc, temp_LMT_oddperc, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMToddperc.vs.CMToddperc", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      
      wilcox.test.temp = wilcox.test(temp_CMT_even, temp_LMT_even, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTeven.vs.CMTeven", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      wilcox.test.temp = wilcox.test(temp_CMT_evenperc, temp_LMT_evenperc, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTevenperc.vs.CMTevenperc", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      
      # CMT vs OMT
      wilcox.test.temp = wilcox.test(temp_CMT, temp_OMT, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMT.vs.CMT", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      
      wilcox.test.temp = wilcox.test(temp_CMT_perc, temp_OMT_perc, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTperc.vs.CMTperc", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      ## by year
      wilcox.test.temp = wilcox.test(temp_CMT_odd, temp_OMT_odd, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTodd.vs.CMTodd", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      wilcox.test.temp = wilcox.test(temp_CMT_oddperc, temp_OMT_oddperc, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMToddperc.vs.CMToddperc", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      
      wilcox.test.temp = wilcox.test(temp_CMT_even, temp_OMT_even, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTeven.vs.CMTeven", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      wilcox.test.temp = wilcox.test(temp_CMT_evenperc, temp_OMT_evenperc, conf.int=TRUE)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTevenperc.vs.CMTevenperc", p = wilcox.test.temp$p.value, estimate = wilcox.test.temp$estimate, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
    } else {
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMT.vs.CMT", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTperc.vs.CMTperc", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTodd.vs.CMTodd", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMToddperc.vs.CMToddperc", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTeven.vs.CMTeven", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "LMTevenperc.vs.CMTevenperc", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMT.vs.CMT", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTperc.vs.CMTperc", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTodd.vs.CMTodd", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMToddperc.vs.CMToddperc", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTeven.vs.CMTeven", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
      treatment_sign = treatment_sign %>% 
        add_row(trait = alltraits[i], depth = depths[j], comparison = "OMTevenperc.vs.CMTevenperc", p = NA, n_scaffoldswtrait = n_scaffoldswtrait, n_scaffoldswtraitfromdepth = n_scaffoldswtraitfromdepth)
    }
    #my_comparisons <- list( c("CMT", "OMT"), c("CMT", "LMT") )
    ## CMT vs OMT and CMT vs LMT
    #p <- ggboxplot(temp, x = "treatment", y = "coverage", 
    #               #color = "treatment", palette = "jco", 
    #               title = paste0("depth: ", depths[j], " (n=", n_scaffoldswtraitfromdepth, " scaffolds)"), 
    #               color = "#00000033", fill = "white", size = 0.5, 
    #               add.params = list(color = "#00000033", size= 0.5),
    #               add = "jitter", ylab = "scaffold coverage", ylim = c(0, 380)) +
    #      font("title", size = 11, color = "red", face = "bold.italic")
#
    #p <- p + stat_compare_means(label = "p", method = "wilcox.test", 
    #                            method.args = list(alternative = "two.sided"), 
    #                            comparisons = my_comparisons, show.legend = F)
    #assign(paste0("p_", sub("-", "", depths[j])), p)
  }
  #p.grid = arrangeGrob(p_015, p_1530, p_3060, p_60100, 
  #                      ncol = 1,
  #                      top = textGrob(alltraits[i], gp = gpar(fontsize=14, fontface = "bold"))
  #                      )
  #alltraits.plotlist[[i]] = p.grid
  #ggsave(p.grid, 
  #       width = 4, height = 16,
  #       filename = file.path(base, "results/", paste0(sub("/", "", alltraits[i]), ".pdf")))
}

#treatment_sign_wide1 = treatment_sign %>%
#  tidyr::pivot_wider(names_from = c("comparison"), values_from = c("p", "estimate"))

treatment_sign_wide = treatment_sign %>%
  tidyr::pivot_wider(names_from = c("comparison"), values_from = c("p", "estimate")) %>%
  #dplyr::mutate_at(c("p_LMT.vs.CMT", "p_LMTperc.vs.CMTperc", "p_OMT.vs.CMT", "p_OMTperc.vs.CMTperc"), p2signiflevels) %>%
  dplyr::mutate(`plevel_LMT.vs.CMT` = unlist(lapply(`p_LMT.vs.CMT`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_LMTperc.vs.CMTperc` = unlist(lapply(`p_LMTperc.vs.CMTperc`, p2signiflevels))) %>% 
  dplyr::mutate(`plevel_LMTodd.vs.CMTodd` = unlist(lapply(`p_LMTodd.vs.CMTodd`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_LMToddperc.vs.CMToddperc` = unlist(lapply(`p_LMToddperc.vs.CMToddperc`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_LMTeven.vs.CMTeven` = unlist(lapply(`p_LMTeven.vs.CMTeven`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_LMTevenperc.vs.CMTevenperc` = unlist(lapply(`p_LMTevenperc.vs.CMTevenperc`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_OMT.vs.CMT` = unlist(lapply(`p_OMT.vs.CMT`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_OMTperc.vs.CMTperc` = unlist(lapply(`p_OMTperc.vs.CMTperc`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_OMTodd.vs.CMTodd` = unlist(lapply(`p_OMTodd.vs.CMTodd`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_OMToddperc.vs.CMToddperc` = unlist(lapply(`p_OMToddperc.vs.CMToddperc`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_OMTeven.vs.CMTeven` = unlist(lapply(`p_OMTeven.vs.CMTeven`, p2signiflevels))) %>%
  dplyr::mutate(`plevel_OMTevenperc.vs.CMTevenperc` = unlist(lapply(`p_OMTevenperc.vs.CMTevenperc`, p2signiflevels))) %>%
  dplyr::mutate(`estimate_LMT.vs.CMT` = as.character(unlist(lapply(`estimate_LMT.vs.CMT`, estimate2sign)))) %>%
  dplyr::mutate(`estimate_LMTperc.vs.CMTperc` = as.character(unlist(lapply(`estimate_LMTperc.vs.CMTperc`, estimate2sign)))) %>% 
  dplyr::mutate(`estimate_LMTodd.vs.CMTodd` = as.character(unlist(lapply(`estimate_LMTodd.vs.CMTodd`, estimate2sign)))) %>%
  dplyr::mutate(`estimate_LMToddperc.vs.CMToddperc` = as.character(unlist(lapply(`estimate_LMToddperc.vs.CMToddperc`, estimate2sign)))) %>% 
  dplyr::mutate(`estimate_LMTeven.vs.CMTeven` = as.character(unlist(lapply(`estimate_LMTeven.vs.CMTeven`, estimate2sign)))) %>%
  dplyr::mutate(`estimate_LMTevenperc.vs.CMTevenperc` = as.character(unlist(lapply(`estimate_LMTevenperc.vs.CMTevenperc`, estimate2sign)))) %>%  
  dplyr::mutate(`estimate_OMT.vs.CMT` = as.character(unlist(lapply(`estimate_OMT.vs.CMT`, estimate2sign)))) %>%
  dplyr::mutate(`estimate_OMTperc.vs.CMTperc` = as.character(unlist(lapply(`estimate_OMTperc.vs.CMTperc`, estimate2sign)))) %>%
  dplyr::mutate(`estimate_OMTodd.vs.CMTodd` = as.character(unlist(lapply(`estimate_OMTodd.vs.CMTodd`, estimate2sign)))) %>%
  dplyr::mutate(`estimate_OMToddperc.vs.CMToddperc` = as.character(unlist(lapply(`estimate_OMToddperc.vs.CMToddperc`, estimate2sign)))) %>% 
  dplyr::mutate(`estimate_OMTeven.vs.CMTeven` = as.character(unlist(lapply(`estimate_OMTeven.vs.CMTeven`, estimate2sign)))) %>%
  dplyr::mutate(`estimate_OMTevenperc.vs.CMTevenperc` = as.character(unlist(lapply(`estimate_OMTevenperc.vs.CMTevenperc`, estimate2sign)))) %>%
  dplyr::select(c("trait", "depth", "n_scaffoldswtrait", "n_scaffoldswtraitfromdepth", matches("(^plevel.*perc)"), matches("(^estimate.*perc)")))   

treatment_sign_wide = treatment_sign_wide %>%
  dplyr::mutate(`select.significant` = case_when(`plevel_LMTperc.vs.CMTperc` %in% c("***", "****") | `plevel_OMTperc.vs.CMTperc` %in% c("***", "****") ~ "1",
                                          TRUE ~ "0"
                                          )) %>%
  dplyr::select(c("trait", "depth", "n_scaffoldswtrait", "n_scaffoldswtraitfromdepth", "select.significant", everything()))

write.table(treatment_sign_wide, file.path(outdir, "microtrait-scaffolds-treatment-enrichment_results.xls"), row.names = F, col.names = T, sep = "\t", quote = FALSE)
saveRDS(treatment_sign_wide, file.path(outdir, "microtrait-scaffolds-treatment-enrichment_results.rds"))



for(i in 1:nrow(treatment_sign_wide)) {
  if(treatment_sign_wide[i,] %>% pull(`select.significant`) == 0) {next}
  cat(i, "\n")}
  row.trait = treatment_sign_wide[i, "trait"] %>% pull
  row.depth = treatment_sign_wide[i, "depth"] %>% pull

  scaffoldswtrait = scaffold2trait %>% 
    dplyr::filter(`trait-display-short` == row.trait) %>% 
    dplyr::select(c("scaffold", "n", "trait-type", "length")) %>% 
    dplyr::mutate(`assembly-depth` = stringr::str_replace(scaffold, "(.*)_(.*)_(.*)_(.*)_(.*)_(.*)_(.*)$", "\\4"))

  scaffoldswtrait = scaffoldswtrait %>% 
    dplyr::filter(length >= 2000)
  n_scaffoldswtrait = scaffoldswtrait %>% dplyr::select(scaffold) %>% n_distinct
  
  scaffoldswtraitfromdepth = scaffoldswtrait %>%
      dplyr::filter(`assembly-depth` == treatment_sign_wide[i, "depth"] %>% pull)
  n_scaffoldswtraitfromdepth = scaffoldswtraitfromdepth %>% dplyr::select(scaffold) %>% n_distinct

  temp = scaffoldswtraitfromdepth %>%
           dplyr::inner_join(scaffold2coverage, by = c("scaffold" = "contig")) %>% 
           dplyr::filter(depth == row.depth) %>%
           dplyr::select(c("coverage", "coverage_perc", "treatment", "year", "n", "trait-type")) %>%
           dplyr::mutate(coverage_weighted = `coverage`*`n`,
                         coverage_perc_weighted = `coverage_perc`*`n`)

  my_comparisons <- list( c("OMT", "CMT"), c("LMT", "CMT") )
  ## CMT vs OMT and CMT vs LMT
  p <- ggboxplot(temp, x = "treatment", y = "coverage_perc_weighted", 
                 #color = "treatment", palette = "jco", 
                 title = paste0("depth: ", row.depth, " (n=", n_scaffoldswtraitfromdepth, " scaffolds)"), 
                 color = "#00000033", fill = "white", size = 0.5, 
                 add.params = list(color = "#00000033", size= 0.5),
                 add = "jitter", ylab = "scaffold coverage") + # ylim = c(0, 100)) +
        font("title", size = 11, color = "red", face = "bold.italic")

  p <- p + stat_compare_means(label = "p", method = "wilcox.test", 
                                method.args = list(alternative = "two.sided"), 
                                comparisons = my_comparisons, show.legend = F)
  ggsave(p.grid, 
         width = 4, height = 16,
         filename = file.path(base, "results/", paste0(sub("/", "", alltraits[i]), ".pdf")))


  #assign(paste0("p_", sub("-", "", depths[j])), p)
}



