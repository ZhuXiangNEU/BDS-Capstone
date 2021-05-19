rm(list=ls())

library(shiny)
library(maplet)
library(tidyverse)
library(DT)
library(plotly)
library(shinyWidgets)



########################################################################## 
########################   Data Exploration  #############################
########################################################################## 

#### Initialize ----
#### Initialize ----
library(maplet)
library(tidyverse)
#### Loading and preprocessing ----
file_data <- system.file("extdata", "example_data/simulated_data.xlsx", package = "maplet")
D <-
    # validate checksum
    mt_load_checksum(file=file_data, checksum = "80afcd72481c6cf3dcf83342e3513699") %>%
    # load data
    mt_load_xls(file=file_data, sheet="data", samples_in_row=T, id_col="sample") %>%
    # load metabolite (rowData) annotations
    mt_anno_xls(file=file_data, sheet="metinfo",anno_type="features", anno_id_col="name", data_id_col = "name") %>%
    # load clinical (colData) annotations
    mt_anno_xls(file=file_data, sheet="clin", anno_type="samples", anno_id_col ="sample", data_id_col ="sample") %>%
    # # log assay dimensions and number of columns for both metabolite and clincial annotations
    mt_reporting_data() %>%
    # generate variables
    mt_anno_mutate(anno_type = "samples", col_name = "outcome1", term = ifelse(Diagnosis==1,rnorm(1,mean=0, sd=1), rnorm(1,mean=0.5, sd=1))) %>%
    mt_anno_mutate(anno_type = "samples", col_name = "outcome2", term = Age+rnorm(1)) %>%
    # heading
    mt_reporting_heading(heading = "Data Clean-up", lvl = 1) %>%
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
    # heading for html file
    mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
    # heading for html file
    mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
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
    # heading for html file
    mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
    # plot sample boxplots
    mt_plots_sample_boxplot(color=Diagnosis, title = "Original", plot_logged = T) %>%
    # apply batch correction
    mt_pre_batch_median(batch_col = "BOX.NUMBER") %>%
    # plot sample boxplots after batch correction
    mt_plots_sample_boxplot(color=Diagnosis, title = "After batch correction", plot_logged = T) %>%
    # normalize abundances using probabilistic quotient
    mt_pre_norm_quot(feat_max = 0.2, ref_samples = Diagnosis==0) %>%
    # show dilution plot
    mt_plots_dilution_factor(in_col="Diagnosis") %>%
    # plot sample boxplots after normalization
    mt_plots_sample_boxplot(color=Diagnosis, title = "After normalization", plot_logged = T) %>%
    # log transform
    mt_pre_trans_log() %>%
    # impute missing values using knn
    mt_pre_impute_knn() %>%
    # plot sample boxplot after imputation
    mt_plots_sample_boxplot(color=Diagnosis, title = "After imputation", plot_logged = T) %>%
    # outlier detection (univariate)
    mt_pre_outlier_detection_univariate() %>%
    # print infos about dataset
    mt_reporting_data() %>%
    # heading for html file
    mt_reporting_heading(heading = "Get Pathway Annotations", lvl = 1) %>%
    # get KEGG ids from HMDB ids
    mt_anno_hmdb_to_kegg(in_col = "HMDb", out_col = "KEGG_ids") %>%
    # get pathway annotations
    #   alternative functions: mt_anno_pathways_xls, mt_anno_pathways_graphite, mt_anno_pathways_uniprot
    mt_anno_pathways_hmdb(in_col = "HMDb", out_col = "pathway", pwdb_name = "KEGG") %>%
    # remove redundant
    mt_anno_pathways_remove_redundant(feat_col = "KEGG_ids", pw_col = "pathway") %>%
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
#### Differential analysis ----
D %<>%
    # heading for html file
    mt_reporting_heading(heading = "Missingness analysis", lvl = 1) %>%
    # compute Fisher's exact test
    mt_stats_univ_missingness(in_col="Diagnosis", stat_name="missingness") %>%
    # create p-value qq plot
    mt_plots_pval_qq(stat_name = "missingness") %>%
    # apply multiple testing correction
    mt_post_multtest(stat_name="missingness", method="BH") %>%
    # heading for html file
    mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
    # heading for html file
    mt_reporting_heading(heading = "Age analysis", lvl = 2) %>%
    # # analysis
    # p_diff_analysis(varname="Age", stat_name= "Age", alpha = 0.05, box_scatter = "scatter") %>%
    # linear model
    mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","Age")), stat_name = "Age") %>%
    # add multiple testing correction
    mt_post_multtest(stat_name = "Age", method = "BH") %>%
    # add stats logging
    mt_reporting_stats(stat_name = "Age", stat_filter = p.adj < 0.05) %>%
    # volcano plot as overview of results
    mt_plots_volcano(stat_name = "Age",
                     x = statistic,
                     feat_filter = p.adj < 0.05,
                     color = p.adj < 0.05) %>%
    # scatter plot
    mt_plots_box_scatter(stat_name = "Age",
                         x = Age,
                         plot_type = "scatter",
                         feat_filter = p.adj < 0.05, 
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
    # barplot
    mt_plots_stats_pathway_bar(stat_list = "Age",
                               feat_filter = p.adj < 0.05,
                               group_col = "SUB_PATHWAY",
                               color_col = "SUPER_PATHWAY") %>%
    # heading for html file
    mt_reporting_heading(heading = "outcome1 analysis", lvl = 2) %>%
    # # analysis
    # p_diff_analysis(varname="outcome1", stat_name= "outcome1", alpha = 0.1, box_scatter = "scatter") %>%
    # linear model
    mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","outcome1")), stat_name = "outcome1") %>%
    # add multiple testing correction
    mt_post_multtest(stat_name = "outcome1", method = "BH") %>%
    # add stats logging
    mt_reporting_stats(stat_name = "outcome1", stat_filter = p.adj < 0.1) %>%
    # volcano plot as overview of results
    mt_plots_volcano(stat_name = "outcome1",
                     x = statistic,
                     feat_filter = p.adj < 0.1,
                     color = p.adj < 0.1) %>%
    # scatter plot
    mt_plots_box_scatter(stat_name = "outcome1",
                         x = outcome1,
                         plot_type = "scatter",
                         feat_filter = p.adj < 0.1, 
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
    # barplot
    mt_plots_stats_pathway_bar(stat_list = "outcome1",
                               feat_filter = p.adj < 0.1,
                               group_col = "SUB_PATHWAY",
                               color_col = "SUPER_PATHWAY") %>%
    # heading for html file
    mt_reporting_heading(heading = "outcome2 analysis", lvl = 2) %>%
    # # analysis
    # p_diff_analysis(varname="outcome2", stat_name= "outcome2", alpha = 0.1, box_scatter = "scatter") %>%
    # linear model
    mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","outcome2")), stat_name = "outcome2") %>%
    # add multiple testing correction
    mt_post_multtest(stat_name = "outcome2", method = "BH") %>%
    # add stats logging
    mt_reporting_stats(stat_name = "outcome2", stat_filter = p.adj < 0.1) %>%
    # volcano plot as overview of results
    mt_plots_volcano(stat_name = "outcome2",
                     x = statistic,
                     feat_filter = p.adj < 0.1,
                     color = p.adj < 0.1) %>%
    # scatter plot
    mt_plots_box_scatter(stat_name = "outcome2",
                         x = outcome2,
                         plot_type = "scatter",
                         feat_filter = p.adj < 0.1, 
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
    # barplot
    mt_plots_stats_pathway_bar(stat_list = "outcome2",
                               feat_filter = p.adj < 0.1,
                               group_col = "SUB_PATHWAY",
                               color_col = "SUPER_PATHWAY") %>%
    # heading for html file
    mt_reporting_heading(heading = "Diagnosis", lvl = 2) %>%
    # # analysis
    # p_diff_analysis(varname="Diagnosis", stat_name= "Diagnosis", alpha = 0.5, box_scatter = "box") %>%
    # linear model
    mt_stats_univ_lm(formula = as.formula(sprintf("~  %s","Diagnosis")), stat_name = "Diagnosis") %>%
    # add multiple testing correction
    mt_post_multtest(stat_name = "Diagnosis", method = "BH") %>%
    # add stats logging
    mt_reporting_stats(stat_name = "Diagnosis", stat_filter = p.adj < 0.5) %>%
    # volcano plot as overview of results
    mt_plots_volcano(stat_name = "Diagnosis",
                     x = statistic,
                     feat_filter = p.adj < 0.5,
                     color = p.adj < 0.5) %>%
    # scatter plot
    mt_plots_box_scatter(stat_name = "Diagnosis",
                         x = Diagnosis,
                         plot_type = "box",
                         feat_filter = p.adj < 0.5, # made very small because otherwise would take an eternity to generate all plots
                         feat_sort = p.value,
                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
    # barplot
    mt_plots_stats_pathway_bar(stat_list = "Diagnosis",
                               feat_filter = p.adj < 0.5,
                               group_col = "SUB_PATHWAY",
                               color_col = "SUPER_PATHWAY") %>%
    # heading for html file
    mt_reporting_heading(heading = "Results Overview", lvl = 1) %>%
    # barplot
    mt_plots_stats_pathway_bar(stat_list = c("Age","outcome1","outcome2","Diagnosis"),
                               feat_filter = p.adj < 0.2,
                               group_col = "SUB_PATHWAY",
                               color_col = "SUPER_PATHWAY") %>%
    {.}

# Extract all the object names for accessor functions
obj_list <- data.frame()
for (i in seq_along(metadata(D)$results)) {
    for (j in seq_along(metadata(D)$results[[i]]$fun)) {
        obj_list[i, j] <- metadata(D)$results[[i]]$fun[j]
    }
}

# Extract all the stat_name
stat_name <- data.frame(stat_name=NA)
for (i in seq_along(metadata(D)$results)) {
    stat_name[i, 1] <- if ('stat_name' %in% names(metadata(D)$results[[i]]$args)) {
        metadata(D)$results[[i]]$args$stat_name
    } else if ('stat_list' %in% names(metadata(D)$results[[i]]$args) & class(metadata(D)$results[[i]]$args$stat_list)!="call") {
        metadata(D)$results[[i]]$args$stat_list
    } else {
        NA
    }
}


# merge object names and stat_name
order_id <- 1:nrow(obj_list)
obj_name <- cbind(order_id, obj_list, stat_name) 
obj_name$stat_name <- ifelse(is.na(obj_name$stat_name), 
                             "(no stat_name)", 
                             obj_name$stat_name)


## SCALE -log10
reverselog_trans <- function (base = exp(1)){
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
}
########################################################################## 
######################## build the Shiny App #############################
########################################################################## 

ui <- fluidPage(
    theme = "bootstrap.css",
    includeCSS("www/style.css"),
    setBackgroundColor("#FFFFFF"),# set canvas background color
    div(style = "padding: 1px 0px; width: '100%'",
        titlePanel(
            title = "",
            windowTitle = "Maplet"
        )
    ),
    # remove shiny "red" warning messages on GUI
    tags$style(type="text/css",
               ".shiny-output-error { visibility: hidden; }",
               ".shiny-output-error:before { visibility: hidden; }"
    ),
    
    navbarPage(
        # embed Maplet logo and title
        title = div(img(src='logo.png',
                        style="float:left; margin-top: -10px; padding-right:10px;padding-bottom:10px",
                        height = 60),
                    "Krumsiek Lab",
                    tags$script(HTML("var header = $('.navbar > .container-fluid');header.append('<div style=\"float:right\"><a href=\"https://github.com/krumsieklab/maplet\"><img src=\"github.png\" alt=\"github\" style=\"float:right;width:33px;height:40px;padding-top:10px;\"> </a></div>');console.log(header)")),
                    br(),
                    tags$script(HTML("var header = $('.navbar > .container-fluid');header.append('<div style=\"float:right\"><a href=\"https://weill.cornell.edu\"><img src=\"WCM.png\" alt=\"logo\" style=\"float:right;height:50px;margin-top: 6px; padding-right:10px; \"> </a></div>');console.log(header)")),
                    windowTitle = "Maplet"),
        # sticky tabs while scrolling main panel
        position = c("fixed-top"), 
        
        # Module 1
        tabPanel("Module 1"),
        
        # Module 2
        tabPanel("Module 2",
                 fluidRow(column(
                     width = 10,
                     selectInput(
                         "mod2.stat",
                         "Select one stat name:",
                         choices = distinct(obj_name[obj_name$V1 == "plots" &
                                                         obj_name$V2 == "stats", ],
                                            stat_name)$stat_name,
                         selected  = "outcome1"
                     ),
                     plotlyOutput("mod2.bar", height = 600),
                     br())), 
                 fluidRow(column(width = 10,
                                 plotlyOutput("mod2.vol", height = 400),
                                 br())),
                 fluidRow(column(width = 10,
                                 plotlyOutput("mod2.box", height = 400)))
                 ), 
        
        
        # Module 3
        tabPanel("Module 3"),
        
        
        # Module 4
        tabPanel("Module 4"),
        
        
        # Module 5
        tabPanel("Module 5"),
        
        
        # Module 6
        tabPanel("Module 6")
        
        
    )
)


server <- function(input, output, session) {
    
    #
    output$mod2.bar <- renderPlotly({
        plots <- mtm_res_get_entries(D, c("plots", "stats"))
        for (i in seq_along(plots)) {
            if (plots[[i]]$args$stat_list == input$mod2.stat) {
                plot <- plots[[i]]$output[[1]]
            }
        }
        ggplotly(plot, source = "sub_bar") %>% layout(dragmode = "lasso")
    })
    
    #
    output$mod2.vol <- renderPlotly({
        d <- event_data("plotly_click", source = "sub_bar")
        if (!is.null(d)) {
            plots <- mtm_res_get_entries(D, c("plots", "stats"))
            for (i in seq_along(plots)) {
                if (plots[[i]]$args$stat_list == input$mod2.stat) {
                    data <- plots[[i]]$output[[1]]$data
                }
            }
            lvls <- rev(levels(data$label))
            label <- lvls[round(as.numeric(d$y))]
            name <- data[data$label == label, ]$name
            
            plots_vol <- mtm_res_get_entries(D, c("plots", "volcano"))
            for (i in seq_along(plots_vol)) {
                if (plots_vol[[i]]$args$stat_name == input$mod2.stat) {
                    data_vol <- plots_vol[[i]]$output[[1]]$data
                }
            }
            
            row_data <- rowData(D) %>% data.frame()
            names <- unlist(row_data[row_data$SUB_PATHWAY == name,]$name)
            data_vol <- data_vol[data_vol$name %in% names, ]
            
            plot <- data_vol %>% ggplot(aes(x = statistic, y = p.value)) +
                geom_point() +
                scale_y_continuous(trans = reverselog_trans(10),
                                   breaks = scales::trans_breaks("log10", function(x) 10^x),
                                   labels = scales::trans_format("log10", scales::math_format(10^.x))) +
                labs(y = "p-value") +
                ggtitle(input$mod2.stat) +
                ggrepel::geom_text_repel(aes(label = name), max.overlaps = Inf)
            ggplotly(plot, source = "sub_vol") %>% layout(dragmode = "lasso")
        }
    })
    
    #
    output$mod2.box <- renderPlotly({
        d1 <- event_data("plotly_click", source = "sub_bar")
        d2 <- event_data("plotly_click", source = "sub_vol")
        if (!is.null(d1) & !is.null(d2)) {
            plots <- mtm_res_get_entries(D, c("plots", "stats"))
            for (i in seq_along(plots)) {
                if (plots[[i]]$args$stat_list == input$mod2.stat) {
                    data <- plots[[i]]$output[[1]]$data
                }
            }
            lvls <- rev(levels(data$label))
            label <- lvls[round(as.numeric(d1$y))]
            name <- data[data$label == label, ]$name
            
            plots_vol <- mtm_res_get_entries(D, c("plots", "volcano"))
            for (i in seq_along(plots_vol)) {
                if (plots_vol[[i]]$args$stat_name == input$mod2.stat) {
                    data_vol <- plots_vol[[i]]$output[[1]]$data
                }
            }
            
            row_data <- rowData(D) %>% data.frame()
            names <- unlist(row_data[row_data$SUB_PATHWAY == name,]$name)
            data_vol <- data_vol[data_vol$name %in% names, ]
            
            
            name2 <- data_vol[as.numeric(d2$pointNumber) + 1, ]$name[1]
            
            plots_box <- mtm_res_get_entries(D, c("plots", "box"))
            for (i in seq_along(plots_box)) {
                if (plots_box[[i]]$args$stat_name == input$mod2.stat) {
                    data_box <- plots_box[[i]]$output[[1]]$data
                }
            }
            data_box <- data_box[data_box$name==name2, ]
            p.value <- signif(mean(data_box$p.value), 3)
            p.adj <- signif(mean(data_box$p.adj), 3)
            plot <- data_box %>%
                ggplot(aes_string(x = input$mod2.stat, y = "value")) +
                geom_point() +
                geom_smooth(method = "lm", se = FALSE) +
                ggtitle(name2) +
                annotate("text",x=-Inf,y=Inf,hjust=0,vjust=2,label=paste0("P-value: ", p.value)) +
                annotate("text",x=-Inf,y=-Inf,hjust=0,vjust=-1,label=paste0("P.adj: ", p.adj))
            ggplotly(plot) %>% layout(dragmode = "lasso")
        }
    })
}


shinyApp(ui = ui, server = server)



