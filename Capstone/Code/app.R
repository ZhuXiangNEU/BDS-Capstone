rm(list=ls())

library(shiny)
library(shinythemes)
library(maplet)
library(tidyverse)
library(DT)
library(dplyr)

######################    Maplet   #########################
#          Jinfeng Lu, Xiang Zhu, Yifan Wu                 #
#                                                          #
############################################################


########################################################################## 
########################   Data Exploration  #############################
########################################################################## 

purrr::zap()

file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "maplet")

# PART 1 - STARTING A METABOTOOLS PIPELINE ----------------------------------------------------

D <-
    # validate checksum
    mt_load_checksum(file=file_data, checksum = "80afcd72481c6cf3dcf83342e3513699") %>%
    # load data - this function loads the assay data only
    #   alternative loading functions: mt_load_metabolon_v1(), mt_load_metabolon_v2(), mt_load_metabolon_lipidomics(),
    #     mt_load_olink(), mt_load_ucd(), mt_load_wcm(), mt_load_nightingale, mt_load_metabolon_new_format()
    mt_load_xls(file=file_data, sheet="data", samples_in_row=T, id_col="sample") %>%
    # load metabolite (rowData) annotations
    mt_anno_xls(file=file_data, sheet="metinfo",anno_type="features", anno_id_col="name", data_id_col = "name") %>%
    # load clinical (colData) annotations
    mt_anno_xls(file=file_data, sheet="clin", anno_type="samples", anno_id_col ="sample", data_id_col ="sample") %>%
    # # log assay dimensions and number of columns for both metabolite and clincial annotations
    mt_reporting_data() %>%
    # start timing
    mt_reporting_tic() %>%
    {.}
# additional functions used at beginning of pipelines:
#   - mt_settings - set global settings for maplet pipeline
#   - mt_load_flag_logged - for flagging a loaded dataset as already log transformed

# PART 2 - DATA CLEANING ----------------------------------------------------

D <- D %>%
    # heading
    mt_reporting_heading(heading = "Data Clean-up", lvl = 1) %>%
    # section text
    mt_reporting_text(text = "Filter samples that are missing values for Diagnosis,add sample annotation column with log10 of
                    PreBioPSA, convert sample annotaiton column Diagnosis to factors,filter metabolites that are missing values
                    for SUB_PATHWAY, log dataset information for this point of the pipeline.") %>%
    # filter samples
    mt_modify_filter_samples(filter = !is.na(Diagnosis)) %>%
    # create additional variable
    mt_anno_mutate(anno_type = "samples", col_name = "PreBioPSALog", term = log10(PreBioPSA)) %>%
    # modify variable to factor
    mt_anno_apply(anno_type = "samples", col_name = "Diagnosis", fun = as.factor) %>%
    # remove metabolites with no pathway annotation
    mt_modify_filter_features(filter = !is.na(SUB_PATHWAY)) %>%
    # log assay dimensions and number of columns for both metabolite and clinical annotations
    mt_reporting_data() %>%
    {.}
# additional data cleaning function:
#   - mt_pre_zero_to_na - for platforms that represent sub-LOD/missing values as zeros


# PART 3.1 - PREPROCESSING: FILTERING MISSING VALUES ----------------------------------------------------

D <- D %>%
    # heading for html file
    mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
    # heading for html file
    mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
    # section text
    mt_reporting_text(text = "Plot percent missingness for each metabolite before filtering, filter out metabolites with >= 50%
                    missingness, plot percent missingness for each metabolite after filtering, add missingness annotation
                    columns to both metabolite and sample annotation data frames.") %>%
    # plot missingness distribution
    mt_plots_missingness(feat_max=0.5) %>%
    # filter metabolites with more than 50% missing values per group
    mt_pre_filter_missingness(feat_max = 0.5, group_col = "Diagnosis") %>%
    # plot missingness distribution after filtering
    mt_plots_missingness(feat_max=0.5) %>%
    # add missingness percentage as annotation to samples (remaining missing)
    mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
    # add missingness percentage as annotation to metabolites
    mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
    {.}


# PART 3.2 - PREPROCESSING: NORMALIZATION ----------------------------------------------------

D <- D %>%
    # heading for html file
    mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
    # section text
    mt_reporting_text(text = "Plot sample boxplots before normalization, apply median batch correction, perform quotient
                    normalization, plot boxplot with dilution factors from quotient normalization, plot sample boxplot after
                    normalization, log transform the data, impute missing data using knn, plot sample boxplot after imputation,
                    detect outliers, log dataset info, write pre-processed data to file.") %>%
    # plot sample boxplots
    mt_plots_sample_boxplot(color=Diagnosis, title = "Original", plot_logged = T) %>%
    # apply batch correction
    #   alternative batch correction function: mt_pre_batch_combat
    mt_pre_batch_median(batch_col = "BOX.NUMBER") %>%
    # plot sample boxplots after batch correction
    mt_plots_sample_boxplot(color=Diagnosis, title = "After batch correction", plot_logged = T) %>%
    # normalize abundances using probabilistic quotient
    #   alternative normalization function: mt_pre_norm_external
    mt_pre_norm_quot(feat_max = 0.2, ref_samples = Diagnosis==0) %>%
    # show dilution plot
    mt_plots_dilution_factor(in_col="Diagnosis") %>%
    # plot sample boxplots after normalization
    mt_plots_sample_boxplot(color=Diagnosis, title = "After normalization", plot_logged = T) %>%
    # log transform
    #   other data transformation functions: mt_pre_trans_exp, mt_pre_trans_relative, mt_pre_trans_scale
    mt_pre_trans_log() %>%
    # impute missing values using knn
    #   alternative imputation functions: mt_pre_impute_min
    mt_pre_impute_knn() %>%
    # plot sample boxplot after imputation
    mt_plots_sample_boxplot(color=Diagnosis, title = "After imputation", plot_logged = T) %>%
    # outlier detection (univariate)
    #   alternative functions: mt_pre_outlier_detection_mahalanobis(), mt_pre_outlier_detection_leverage()
    #   related function: mt_pre_outlier_to_na()
    mt_pre_outlier_detection_univariate() %>%
    # print infos about dataset
    mt_reporting_data() %>%
    # write preprocessed data to Excel file
    #   other writing functions: mt_write_se_rds (save SummarizedExerpiment object)
    mt_write_se_xls(file = "PreprocessedData.xlsx") %>%
    {.}

# Additional pre-processing functions
#   - mt_pre_confounding_correction() - function for correcting confounding variables
#   - mt_pre_confounding_correction_stepwise_aic() - an alterenative function for correcting confounders that uses stepwise aic
# NOTES ON BEST PRACTICES: If incorporated in this pipeline, the above functions would correct for the variable age such that
#     none of the following functions have to take care of those confounders anymore. If this function is included, confounders
#     should not be included in any of the following functions. It is generally agreed that including the confounders in the
#     linear models themselves is preferable to pre-correction.

# PART 4 - GET PATHWAY ANNOTATIONS ----------------------------------------------------

D <- D %>%
    # heading for html file
    mt_reporting_heading(heading = "Get Pathway Annotations", lvl = 1) %>%
    # get KEGG ids from HMDB ids
    mt_anno_hmdb_to_kegg(in_col = "HMDb", out_col = "KEGG_ids") %>%
    # get pathway annotations
    #   alternative functions: mt_anno_pathways_xls, mt_anno_pathways_graphite, mt_anno_pathways_uniprot
    mt_anno_pathways_hmdb(in_col = "HMDb", out_col = "pathway", pwdb_name = "KEGG") %>%
    # remove redundant
    mt_anno_pathways_remove_redundant(feat_col = "KEGG_ids", pw_col = "pathway") %>%
    # write pathway annotations
    mt_write_pathways(file="ExamplePipeline_PathwayAnnotations.xlsx", pw_col = "pathway") %>%
    {.}


# PART 5 - GLOBAL STATISTICS ----------------------------------------------------

D <- D %>%
    # heading for html file
    mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
    # plot PCA
    mt_plots_pca(scale_data = T, title = "scaled PCA - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
    # plot UMAP
    mt_plots_umap(scale_data = T, title = "scaled UMAP - Diagnosis", color=Diagnosis, size=2.5, ggadd=scale_size_identity()) %>%
    # plot heatmap
    mt_plots_heatmap(scale_data = T, annotation_col = c("Diagnosis"), annotation_row = c("SUPER_PATHWAY"),
                     clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
    {.}

# PART 6.1 - STATISTICAL ANALYSIS, OUTCOME: DIAGNOSIS, METHOD: MISSINGNESS ANALYSIS ---------------------------------------

#create another SE object for first analysis branch (missingness & metabolites)
D1 <- D

D1 <- D1 %>%
    # heading for html file
    mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
    # section text
    mt_reporting_text(text = "Perform missingness analysis to determine if NAs significantly accumulate in one of the Diagnosis
                    groups. Adjust output of test using multiple testing correction.") %>%
    # compute Fisher's exact test
    mt_stats_univ_missingness(in_col="Diagnosis", stat_name="missingness") %>%
    # create p-value qq plot
    mt_plots_pval_qq(stat_name = "missingness") %>%
    # apply multiple testing correction
    #   alternative function: mt_post_multtest_effdim
    mt_post_multtest(stat_name="missingness", method="BH") %>%
    {.}


########################################################################## 
########################       Functions     #############################
########################################################################## 

# Extract all the object names for accessor
obj_list <- data.frame()
for (i in seq_along(metadata(D1)$results)) {
    for (j in seq_along(metadata(D1)$results[[i]]$fun)) {
        obj_list[i, j] <- metadata(D1)$results[[i]]$fun[j]
    }
}

# Extract all the stat_name
stat_name <- data.frame(name=NA)
for (i in seq_along(metadata(D1)$results)) {
    if ("stat_name" %in% names(metadata(D1)$results[[i]]$args)){
        stat_name[i, 1] <- metadata(D1)$results[[i]]$args$stat_name
    } else {
        stat_name[i, 1] <- "(no stat_name)"
    }
}
# distinct the not-null values
# stat_name <- distinct(subset(stat_name, !is.na(name)), name)

# merge object names and stat_name
obj_name <- cbind(obj_list, stat_name)

##
get_plots_by_tabs <- function(tab1, tab2) {
    results <- D1 %>% mtm_res_get_entries(c(tab1, tab2))
    
    plots <- list()
    for (i in 1:length(results)) {
        plots[i] <- results[[i]]$output
    }

    return(plots)
}

########################################################################## 
######################## build the Shiny App #############################
########################################################################## 


ui <- shinyUI(
    fluidPage(
        theme = "style.css",
        div(style = "padding: 1px 0px; width: '100%'",
            titlePanel(
                title = "",
                windowTitle = "Maplet"
            )
        ),
        
        navbarPage(
            title = div(img(src='logo.png',
                            style="margin-top: -14px; padding-right:10px;padding-bottom:10px",
                            height = 60),
                        "BDS Capstone Maplet",
                        windowTitle = "Maplet"),
            
            theme = shinytheme("lumen"), 
            
            # Module 1 
            tabPanel(
                "Module 1",
                
                # One tab for each plot/table.
                tabsetPanel(type = "tabs",
                            
                            # Circle-packing plot of ethnicity and gender.
                            tabPanel(
                                "Plots",
                                # Sidebar panel for controls.
                                sidebarPanel(
                                    selectInput(
                                        "mod1.0",
                                        label = "Select stat name:",
                                        choices = c(unique(obj_name$name)),
                                        selected = "(no stat_name)"
                                    ),
                                    uiOutput("mod1.1"),
                                    uiOutput("mod1.2"),
                                    tags$p(
                                        HTML("<b>Hint:</b> Module 1 requires extracting all the result objects."
                                        )),
                                    tags$p(
                                        HTML("Users can assess results in a drop-down menu that offers a list of a statname and a plot type (e.g. “missingness”, “pval”)."
                                        )),
                                    br(),
                                    br(),
                                    # delay the output
                                    actionButton("mod1.go.p", "Update"),
                                    width = 3
                                ),
                                # Main panel with plot.
                                mainPanel(plotOutput("mod1.p"),
                                          width = 9)
                            ),
                            tabPanel("Tables",
                                     
                                     # Sidebar panel for controls.
                                     sidebarPanel(selectInput(
                                         "mod1.3",
                                         label = "Select stat name:",
                                         choices = c(unique(obj_name$name)),
                                         selected = "(no stat_name)"),
                                         br(),
                                         br(),
                                         # delay the output
                                         actionButton("mod1.go.t", "Update"),
                                         width = 3), 
                                     # Main panel with plot.
                                     mainPanel(fluidRow(dataTableOutput("mod1.t")))))),
            
            # Module 2
            tabPanel("Module 2"),
            
            
            # Module 3
            tabPanel("Module 3"),
            
            
            # Module 4
            tabPanel("Module 4"),
            
            
            # Module 5
            tabPanel("Module 5"),
            
            
            # Module 6
            tabPanel("Module 6")
                            
                            
)))


server <- shinyServer(function(input, output) {
    
    output$mod1.1 = renderUI({
        selectInput(
            "mod1.1",
            label = "Select plot type:",
            choices = if (input$"mod1.0" == "(no stat_name)") {
                c(unique(obj_name$V2[obj_name$V1 == "plots"]))
            } else {
                c(unique(obj_name$V2[obj_name$V1 == "plots" & obj_name$name == input$"mod1.0"]))
            },
            selected = ""
        )
    })
    
    output$mod1.2 = renderUI({
        selectInput(
            "mod1.2",
            label = "Select plot number:",
            choices = c(1:length(get_plots_by_tabs("plots", input$"mod1.1"))),
            selected = ""
        )
    })
    
    # reactive expression
    mode1_plot_reactive <- eventReactive(input$mod1.go.p, {
        c(input$mod1.0, input$mod1.1, input$mod1.2)
    })
    
    output$mod1.p <- renderPlot({
        inputs <- mode1_plot_reactive()
        plots <- get_plots_by_tabs("plots", inputs[2])
        plots[as.numeric(inputs[3])]
    })
    
    # reactive expression
    mode1_table_reactive <- eventReactive(input$mod1.go.t, {
        c(input$mod1.3)
    })
    
    output$mod1.t <- renderDataTable({
        inputs <- mode1_table_reactive()
        mtm_get_stat_by_name(D1, inputs[1])
    })
})



shinyApp(ui = ui, server = server)
