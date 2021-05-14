rm(list=ls())

library(shiny)
library(shinyWidgets)
library(maplet)
library(tidyverse)
library(DT)
library(plotly)

######################    Maplet   #########################
#          Jinfeng Lu, Xiang Zhu, Yifan Wu                 #
#                                                          #
############################################################


########################################################################## 
########################   Data Exploration  #############################
########################################################################## 

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
#### Analysis Function ----
p_diff_analysis <- function(D, varname, stat_name, alpha=0.05, box_scatter="box"){
    D %>%
        # linear model
        mt_stats_univ_lm(formula = as.formula(sprintf("~  %s",varname)), stat_name = stat_name) %>%
        # add multiple testing correction
        mt_post_multtest(stat_name = stat_name, method = "BH") %>%
        # add stats logging
        mt_reporting_stats(stat_name = stat_name, stat_filter = p.adj < alpha) %>%
        # volcano plot as overview of results
        mt_plots_volcano(stat_name = stat_name,
                         x = statistic,
                         feat_filter = p.adj < alpha,
                         color = p.adj < alpha) %>%
        # scatter plot
        mt_plots_box_scatter(stat_name = stat_name,
                             x = !!sym(varname),
                             plot_type = box_scatter,
                             feat_filter = p.adj < alpha, # made very small because otherwise would take an eternity to generate all plots
                             feat_sort = p.value,
                             annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
        # barplot
        mt_plots_stats_pathway_bar(stat_list = stat_name,
                                   feat_filter = p.adj < alpha,
                                   group_col = "SUB_PATHWAY",
                                   color_col = "SUPER_PATHWAY") %>%
        {.}
}
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
    # analysis
    p_diff_analysis(varname="Age", stat_name= "Age", alpha = 0.05, box_scatter = "scatter") %>%
    # heading for html file
    mt_reporting_heading(heading = "outcome1 analysis", lvl = 2) %>%
    # analysis
    p_diff_analysis(varname="outcome1", stat_name= "outcome1", alpha = 0.1, box_scatter = "scatter") %>%
    # heading for html file
    mt_reporting_heading(heading = "outcome2 analysis", lvl = 2) %>%
    # analysis
    p_diff_analysis(varname="outcome2", stat_name= "outcome2", alpha = 0.1, box_scatter = "scatter") %>%
    # heading for html file
    mt_reporting_heading(heading = "Diagnosis", lvl = 2) %>%
    # analysis
    p_diff_analysis(varname="Diagnosis", stat_name= "Diagnosis", alpha = 0.5, box_scatter = "box") %>%
    # heading for html file
    mt_reporting_heading(heading = "Results Overview", lvl = 1) %>%
    # barplot
    mt_plots_stats_pathway_bar(stat_list = c("Age","outcome1","outcome2","Diagnosis"),
                               feat_filter = p.adj < 0.2,
                               group_col = "SUB_PATHWAY",
                               color_col = "SUPER_PATHWAY") %>%
    {.}


########################################################################## 
########################       Functions     #############################
########################################################################## 

# Extract all the object names for accessor
obj_list <- data.frame()
for (i in seq_along(metadata(D)$results)) {
    for (j in seq_along(metadata(D)$results[[i]]$fun)) {
        obj_list[i, j] <- metadata(D)$results[[i]]$fun[j]
    }
}

# Extract all the stat_name
stat_name <- c()
for (i in seq_along(metadata(D)$results)) {
    fun <- metadata(D)$results[[i]]
    stat_name[i] <- if ("stat_name" %in% names(fun$args)) {
        if (fun$fun[1] == "plots") {
            if (fun$args$stat_name == "stat_name") {
                col_name <- colnames(fun$output[[1]]$data)[2]
                ifelse(col_name == "term",
                       unique(fun$output[[1]]$data$term),
                       col_name)
            } else {
                fun$args$stat_name
            }
        } else {
            ifelse(is.null(fun$output),
                   "NULL",
                   fun$output$outcome)
        }
    } else {
        "(no stat_name)"
    }
}
# distinct the not-null values
# stat_name <- distinct(subset(stat_name, !is.na(name)), name)

# merge object names and stat_name
#obj_name <- cbind(obj_list, stat_name)
obj_list$name <- stat_name
obj_list <- obj_list[obj_list$name != "NULL",]


##
get_V2_by_tabs <- function(tab1, tab2, tab3 = NULL) {
    
    if (is.null(tab3)) {
        data <- obj_list[obj_list$V1 == tab1 & obj_list$name == tab2,]
        col_name <- unique(data$V2)
    } else {
        data <- obj_list[obj_list$V1 == tab1 &
                             obj_list$name == tab2 &
                             obj_list$V2 == tab3,]
        col_name <- unique(data$V3)
    }

    return(col_name)
}

get_V3_by_tabs <- function(tab1, tab2, tab3) {
    
    data <- obj_list[obj_list$V1 == tab1 & obj_list$name ==tab2, ]
    col_name <- unique(data$V2)
    
    return(col_name)
}


## pca/umap #, FUN = mt_plots_pca
mod3_plots_pca <- function(num_assay = 1,
                           if_scale = FALSE,
                           if_factor = FALSE,
                           col_name) {
    
    X <- t(assay(D, num_assay))
    pca <- stats::prcomp(x = as.matrix(X), center = F, scale = if_scale)
    df <- data.frame(PCA1 = pca$x[, 1], PCA2 = pca$x[, 2], colData(D))
    expvar <- (pca$sdev)^2 / sum(pca$sdev^2)
    
    if (if_factor) {
        df[, col_name] <- factor(df[, col_name])
    }
    
    df %>% ggplot(
           aes_string(
               x = "PCA1",
               y = "PCA2",
               color = col_name
           )) +
        geom_point() +
        xlab(sprintf("PCA1 (%.1f%%)", expvar[1] * 100)) +
        ylab(sprintf("PCA2 (%.1f%%)", expvar[2] * 100)) +
        ggtitle(paste0("PCA ", "color by ", col_name,
                       " (scaled:", if_scale, ", array:", num_assay, ")")) +
        labs(color = col_name)
}

mod3_plots_umap <- function(num_assay = 1,
                            if_scale = FALSE,
                            if_factor = FALSE,
                            col_name) {
    
    X <- t(assay(D, num_assay))
    
    if (if_scale) {
        X <- scale(X)
    }
    
    umap <- umap::umap(d = as.matrix(X), n_neighbors = 15)
    df <- data.frame(comp1 = umap$layout[, 1], comp2 = umap$layout[, 2], colData(D))
    
    if (if_factor) {
        df[, col_name] <- factor(df[, col_name])
    }
    
    df %>% ggplot(
        aes_string(
            x = "comp1",
            y = "comp2",
            color = col_name
        )) +
        geom_point() +
        ggtitle(paste0("Umap ", "color by ", col_name,
                       " (scaled:", if_scale, ", array:", num_assay, ")")) +
        labs(color = col_name)
}

########################################################################## 
######################## build the Shiny App #############################
########################################################################## 

ui <- shinyUI(
    fluidPage(
        theme = "bootstrap.css",
        includeCSS("www/style.css"),
        setBackgroundColor("#FFFFFF"),# set canvas background color
        div(style = "padding: 1px 0px; width: '100%'",
            titlePanel(
                title = "",
                windowTitle = "Maplet"
            )
        ),
        
        navbarPage(
            title = div(img(src='logo.png',
                            style="float:left; margin-top: -14px; padding-right:10px;padding-bottom:10px",
                            height = 60),
                        "BDS Capstone Maplet",
                        tags$script(HTML("var header = $('.navbar > .container-fluid');header.append('<div style=\"float:right\"><a href=\"https://github.com/krumsieklab/maplet\"><img src=\"github.png\" alt=\"github\" style=\"float:right;width:33px;height:40px;padding-top:10px;\"> </a></div>');console.log(header)")),
                        br(),
                        tags$script(HTML("var header = $('.navbar > .container-fluid');header.append('<div style=\"float:right\"><a href=\"https://weill.cornell.edu\"><img src=\"WCM.png\" alt=\"logo\" style=\"float:right;height:50px;padding-top:-14px;\"> </a></div>');console.log(header)")),
                        windowTitle = "Maplet"),
            
            # Module 1 
            tabPanel(
                "Module 1",
                
                # Sidebar panel for controls.
                sidebarPanel(
                    selectInput(
                        "mod1.1",
                        label = "Select stat name:",
                        choices = c(unique(obj_list$name)),
                        selected = "(no stat_name)"
                    ),
                    radioButtons(
                        "mod1.2",
                        label = "Select output type:",
                        choices = list("Plot" = "plots",
                                       "Table" = "stats"),
                        selected = "plots"
                    ),
                    uiOutput("mod1.3"),
                    uiOutput("mod1.4"),
                    tags$p(
                        HTML(
                            "<b>Hint:</b> Module 1 requires extracting all the result objects."
                        )
                    ),
                    tags$p(
                        HTML(
                            "Users can assess results in a drop-down menu that offers a list of a statname and a plot type (e.g. “missingness”, “pval”)."
                        )
                    ),
                    br(),
                    br(),
                    # delay the output
                    actionButton("mod1.go", "Update"),
                    width = 3
                ),
                # Main panel with plot.
                mainPanel(uiOutput("mod1.out"),
                          width = 9)
            ), 
            
            # Module 2
            tabPanel("Module 2"),
            
            
            # Module 3
            tabPanel("Module 3",
                     # Sidebar panel for controls.
                     sidebarPanel(
                         # check box for scale
                         checkboxInput("mod3.scale",
                                       "Scale",
                                       value = FALSE),
                         # check box for factor
                         checkboxInput("mod3.factor",
                                       "Factor",
                                       value = FALSE),
                         # dropdown menu for different assay
                         selectInput("mod3.assay",
                                     "Select one assay:", 
                                     choices = c(1:length(assays(D))),
                                     selected = ""
                         ),
                         # select one plot type
                         radioButtons("mod3.1",
                                     "Select plot type:",
                                     choices = list("PCA" = "pca", 
                                                    "UMAP" = "umap")), 
                         # select one colData column for coloring
                         selectInput(
                             "mod3.2",
                             "Select colData column:",
                             choices = names(colData(D)),
                             selected = ""
                         ), 
                         tags$p(
                             HTML("<b>Hint:</b> Module 3 requires generating an interactive 2D projection of PCA/UMAP."
                             )),
                         tags$p(
                             HTML("It displays a drop-down menu of all colData columns for coloring."
                             )),
                         br(),
                         br(),
                         # delay the output
                         actionButton("mod3.go", "Update"),
                         width = 3), 
                     # Main panel with plot.
                     mainPanel(plotlyOutput("mod3.out"),
                               width = 9)),
            
            
            # Module 4
            tabPanel("Module 4"),
            
            
            # Module 5
            tabPanel("Module 5"),
            
            
            # Module 6
            tabPanel("Module 6")
                            
                            
)))


server <- shinyServer(function(input, output) {
    
    output$mod1.3 = renderUI({
        selectInput(
            "mod1.3",
            label = "Select the object:",
            choices = get_V2_by_tabs(input$mod1.2, input$mod1.1),
            selected = ""
        )
    })
    
    output$mod1.4 = renderUI({
        if (input$mod1.2 == "plots") {
            return()
        } else {
            selectInput(
                "mod1.4",
                label = "Select the object2:",
                choices = get_V2_by_tabs(input$mod1.2,#stats
                                         input$mod1.1,#Diagnosis
                                         input$mod1.3),#univ
                selected = ""
            )
        }
    })
    
    
    # reactive expression
    mod1_reactive <- eventReactive(input$mod1.go, {
        c(input$mod1.1, input$mod1.2, input$mod1.3, input$mod1.4)
    })
    
    output$mod1.p <- renderUI({
        inputs <- mod1_reactive()
        plots <- D %>% mtm_res_get_entries(c("plots", inputs[3]))
        
        # there are multiple plots
        len_i <- length(plots)
        # some plots have multiple objects
        len_j <- length(plots[[1]]$output)
        # name every plot object in UI
        mod1_plot_output_list <- lapply(1:(len_i * len_j), function(i) {
            plotname <- paste("Plot", i, sep = "")
            plotOutput(plotname)
        })
        
        # Convert the list to a tagList - this is necessary for the list of items
        # to display properly.
        do.call(tagList, mod1_plot_output_list)
    })
    
    # Call renderPlot for each one. Plots are only actually generated when they
    # are visible on the web page.
    for (i in 1:10) {
        # Need local so that each item gets its own number. Without it, the value
        # of i in the renderPlot() will be the same across all instances, because
        # of when the expression is evaluated.
        local({
            my_i <- i
            plotname <- paste("Plot", my_i, sep="")
            
            # render plots
            output[[plotname]] <- renderPlot({
                inputs <- mod1_reactive()
                plots <- D %>% mtm_res_get_entries(c("plots", inputs[3]))
                # there are multiple plots
                len_i <- length(plots)
                # some plots have multiple objects
                len_j <- length(plots[[1]]$output)
                # locate the row in the `plots`
                row_n <- ceiling(my_i/len_i)
                # locate the column in the `plots`
                col_n <- ifelse((my_i %% len_i)==0, len_j, (my_i %% len_i))
                
                plots[[row_n]]$output[col_n]
            })
        })
    }
    
    output$mod1.t <- renderDataTable({
        inputs <- mod1_reactive()
        tb <- D %>% mtm_res_get_entries(c(inputs[2], inputs[3], inputs[4]))
        
        tb[[1]]$output$table
    })
    
    output$mod1.out <- renderUI({
        switch(
            input$mod1.2,
            "plots" = uiOutput("mod1.p"),
            "stats" = fluidRow(dataTableOutput("mod1.t"))
        )
    })
    
    # reactive expression
    mod3_reactive <- eventReactive(input$mod3.go, {
        c(input$mod3.scale,
          input$mod3.factor,
          input$mod3.assay,
          input$mod3.1,
          input$mod3.2)
    })
    
    output$mod3.out <- renderPlotly({
        inputs <- mod3_reactive()
        plot <-
            switch(inputs[4],
                   "pca" = ggplotly(
                       mod3_plots_pca(
                           as.numeric(inputs[3]),
                           as.logical(inputs[1]),
                           as.logical(inputs[2]),
                           inputs[5]
                           
                       )
                   ),
                   "umap" = ggplotly(
                       mod3_plots_umap(
                           as.numeric(inputs[3]),
                           as.logical(inputs[1]),
                           as.logical(inputs[2]),
                           inputs[5]
                       )
                   ))
        
        plot %>%
            layout(
                legend = list(
                    orientation = "h",
                    # show entries horizontally
                    xanchor = "center",
                    # use center of legend as anchor
                    x = .5,
                    y = -.2
                ),
                autosize = TRUE
            )
    })
    
})


shinyApp(ui = ui, server = server)



