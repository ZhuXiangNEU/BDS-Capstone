#### Initialize ----
library(maplet)
library(tidyverse)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#### Loading and preprocessing ----
file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "maplet")
D <-
  mt_load_checksum(file=file_data, checksum = "80afcd72481c6cf3dcf83342e3513699") %>%
  mt_load_xls(file=file_data, sheet="data", samples_in_row=T, id_col="sample") %>%
  mt_anno_xls(file=file_data, sheet="metinfo",anno_type="features", anno_id_col="name", data_id_col = "name") %>%
  mt_anno_xls(file=file_data, sheet="clin", anno_type="samples", anno_id_col ="sample", data_id_col ="sample") %>%
  mt_reporting_data() %>%
  mt_anno_mutate(anno_type = "samples", col_name = "outcome1", term = ifelse(Diagnosis==1,rnorm(1,mean=0, sd=1), rnorm(1,mean=0.5, sd=1))) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "outcome2", term = Age+rnorm(1)) %>%
  {.}

set.seed(1)
D <- D[sample(x=1:nrow(D),size = 300,replace = F),
       sample(x=1:ncol(D),size = 300,replace = F)] %>%
  mt_reporting_heading(heading = "Data Clean-up", lvl = 1) %>%
  mt_modify_filter_samples(filter = !is.na(Diagnosis)) %>%
  mt_anno_mutate(anno_type = "samples", col_name = "PreBioPSALog", term = log10(PreBioPSA)) %>%
  mt_anno_apply(anno_type = "samples", col_name = "Diagnosis", fun = as.factor) %>%
  mt_modify_filter_features(filter = !is.na(SUB_PATHWAY)) %>%
  mt_reporting_data() %>%
  mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
  mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
  mt_plots_missingness(feat_max=0.5) %>%
  mt_pre_filter_missingness(feat_max = 0.5, group_col = "Diagnosis") %>%
  mt_plots_missingness(feat_max=0.5) %>%
  mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
  mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
  mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
  mt_plots_sample_boxplot(color=Diagnosis, title = "Original", plot_logged = T) %>%
  mt_pre_batch_median(batch_col = "BOX.NUMBER") %>%
  mt_plots_sample_boxplot(color=Diagnosis, title = "After batch correction", plot_logged = T) %>%
  mt_pre_norm_quot(feat_max = 0.2, ref_samples = Diagnosis==0) %>%
  mt_plots_dilution_factor(in_col="Diagnosis") %>%
  mt_plots_sample_boxplot(color=Diagnosis, title = "After normalization", plot_logged = T) %>%
  mt_pre_trans_log() %>%
  mt_pre_impute_knn() %>%
  mt_plots_sample_boxplot(color=Diagnosis, title = "After imputation", plot_logged = T) %>%
  mt_pre_outlier_detection_univariate() %>%
  mt_reporting_data() %>%
  mt_reporting_heading(heading = "Get Pathway Annotations", lvl = 1) %>%
  mt_anno_hmdb_to_kegg(in_col = "HMDb", out_col = "KEGG_ids") %>%
  mt_anno_pathways_hmdb(in_col = "HMDb", out_col = "pathway", pwdb_name = "KEGG") %>%
  mt_anno_pathways_remove_redundant(feat_col = "KEGG_ids", pw_col = "pathway") %>%
  mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
  mt_plots_pca(scale_data = T, title = "scaled PCA - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_umap(scale_data = T, title = "scaled UMAP - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
  mt_plots_heatmap(scale_data = T, annotation_col = c("Diagnosis"), annotation_row = c("SUPER_PATHWAY"),
                   clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
  {.}

#### Differential analysis ----

D %<>%
  mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
  mt_stats_univ_missingness(in_col="Diagnosis", stat_name="Fisher's test - Missingness") %>%
  mt_plots_pval_qq(stat_name = "Fisher's test - Missingness") %>%
  mt_post_multtest(stat_name="Fisher's test - Missingness", method="BH") %>%
  mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
  mt_reporting_heading(heading = "Age analysis", lvl = 2) %>%
  mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","Age")), stat_name = "Age Regression") %>%
  mt_post_multtest(stat_name = "Age Regression", method = "BH") %>%
  mt_reporting_stats(stat_name = "Age Regression", stat_filter = p.adj < 0.1) %>%
  mt_plots_volcano(stat_name = "Age Regression",
                   x = statistic,
                   feat_filter = p.adj < 0.1,
                   color = p.adj < 0.1) %>%
  mt_plots_box_scatter(stat_name = "Age Regression",
                       x = Age,
                       plot_type = "scatter",
                       feat_filter = p.adj < 0.1, 
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
  mt_plots_stats_pathway_bar(stat_list = "Age Regression",
                             y_scale = "count",
                             feat_filter = p.adj < 0.1,
                             group_col = "SUB_PATHWAY",
                             color_col = "SUPER_PATHWAY") %>%
  mt_reporting_heading(heading = "outcome1 analysis", lvl = 2) %>%
  mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","outcome1")), stat_name = "Comparison Outcome1") %>%
  mt_post_multtest(stat_name = "Comparison Outcome1", method = "BH") %>%
  mt_reporting_stats(stat_name = "Comparison Outcome1", stat_filter = p.adj < 0.2) %>%
  mt_plots_volcano(stat_name = "Comparison Outcome1",
                   x = statistic,
                   feat_filter = p.adj < 0.2,
                   color = p.adj < 0.2) %>%
  mt_plots_box_scatter(stat_name = "Comparison Outcome1",
                       x = outcome1,
                       plot_type = "box",
                       feat_filter = p.adj < 0.2, 
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
  mt_plots_stats_pathway_bar(stat_list = "Comparison Outcome1",
                             y_scale = "count",
                             feat_filter = p.adj < 0.2,
                             group_col = "SUB_PATHWAY",
                             color_col = "SUPER_PATHWAY") %>%
  mt_reporting_heading(heading = "outcome2 analysis", lvl = 2) %>%
  mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","outcome2")), stat_name = "Comparison Outcome2") %>%
  mt_post_multtest(stat_name = "Comparison Outcome2", method = "BH") %>%
  mt_reporting_stats(stat_name = "Comparison Outcome2", stat_filter = p.adj < 5E-4) %>%
  mt_plots_volcano(stat_name = "Comparison Outcome2",
                   x = statistic,
                   feat_filter = p.adj < 5E-4,
                   color = p.adj < 5E-4) %>%
  mt_plots_box_scatter(stat_name = "Comparison Outcome2",
                       x = outcome2,
                       plot_type = "scatter",
                       feat_filter = p.adj < 5E-4, 
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
  mt_plots_stats_pathway_bar(stat_list = "Comparison Outcome2",
                             y_scale = "count",
                             feat_filter = p.adj < 5E-4,
                             group_col = "SUB_PATHWAY",
                             color_col = "SUPER_PATHWAY") %>%
  mt_reporting_heading(heading = "Diagnosis", lvl = 2) %>%
  mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","Diagnosis")), stat_name = "Diagnosis Analysis") %>%
  mt_post_multtest(stat_name = "Diagnosis Analysis", method = "BH") %>%
  mt_reporting_stats(stat_name = "Diagnosis Analysis", stat_filter = p.adj < 0.3) %>%
  mt_plots_volcano(stat_name = "Diagnosis Analysis",
                   x = statistic,
                   feat_filter = p.adj < 0.3,
                   color = p.adj < 0.3) %>%
  mt_plots_box_scatter(stat_name = "Diagnosis Analysis",
                       x = Diagnosis,
                       plot_type = "box",
                       feat_filter = p.adj < 0.3, # made very small because otherwise would take an eternity to generate all plots
                       feat_sort = p.value,
                       annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
  mt_plots_stats_pathway_bar(stat_list = "Diagnosis Analysis",
                             y_scale = "count",
                             feat_filter = p.adj < 0.3,
                             group_col = "SUB_PATHWAY",
                             color_col = "SUPER_PATHWAY") %>%
  
  mt_reporting_heading(heading = "Results Overview", lvl = 1) %>%
  mt_plots_stats_pathway_bar(stat_list = c("Age Regression","Comparison Outcome1",
                                           "Comparison Outcome2","Diagnosis Analysis"),
                             y_scale = "count",
                             feat_filter = p.adj < 0.05,
                             group_col = "SUB_PATHWAY",
                             color_col = "SUPER_PATHWAY") %>%
  {.}