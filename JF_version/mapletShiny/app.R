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

# define PCA output function for mod3 referring 'mt_plots_pca'
mod3_plots_pca <- function(D, title = "PCA", 
                           ## scale argument
                           scale_data,
                           ## color argument
                           color,
                           categorizing,
                           pc1 = 1, pc2 = 2, 
                           data_type,
                           ...){
  X = t(assay(D))
  if (any(is.na(X))) 
    stop("Data matrix for PCA cannot contain NAs")
  if (!(data_type %in% c("scores", "loadings"))) 
    stop("Show must be either 'scores' or 'loadings'")
  # scale scores if scale checkbox=T
  if (scale_data) 
    X <- scale(X)
  pca <- stats::prcomp(x = as.matrix(X), center = F, scale = F)
  expvar <- (pca$sdev)^2/sum(pca$sdev^2)
  ## reactivate axis labels
  pc1name <- sprintf("PC%d (%.1f%%)", pc1, expvar[pc1] * 
                       100)
  pc2name <- sprintf("PC%d (%.1f%%)", pc2, expvar[pc2] * 
                       100)
  
  if (data_type == "scores") {
    df <- data.frame(x = pca$x[, pc1], 
                    y = pca$x[, pc2], 
                    colData(D)
                    )
    colnames(df)[1:2] <- c(sprintf("PC%d", pc1), 
                           sprintf("PC%d", pc2)
                           )
    ## reactivate plot title
    plot_title <- paste0(ifelse(scale_data,
                                "Scaled ", 
                                "Non-scaled "),
                         title, " - ", color)
    ## categorize coloring if color checkbox=T
    if(categorizing){
      df[, color] <- factor(df[, color])
    }
    # draw ggplot
    p <- ggplot(data = df, 
                aes_string(
                  x = sprintf("PC%d", pc1), 
                  y = sprintf("PC%d", pc2),
                  color = color)) + 
      geom_point() + 
      xlab(pc1name) + 
      ylab(pc2name) + 
      ggtitle(plot_title) +
      ## reactive legend title
      labs(color = color)
  } else {
    df = data.frame(x = pca$rotation[, pc1], 
                    y = pca$rotation[, pc2], 
                    rowData(D)
    )
    colnames(df)[1:2] <- c(sprintf("PC%d", pc1), 
                           sprintf("PC%d", pc2)
    )
    ## reactivate plot title
    plot_title <- paste0(ifelse(scale_data, 
                                "Scaled ", 
                                "Non-scaled "),
                         title)
    # draw ggplot
    p <- ggplot(data = df, 
                aes_string(
                  x = sprintf("PC%d", pc1), 
                  y = sprintf("PC%d", pc2)
                )) + 
      geom_point() + 
      xlab(pc1name) + 
      ylab(pc2name) + 
      ggtitle(plot_title)
  }
  
  # draw plotly
  ggplotly(p) %>%
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = .5, ## set position of legend
                         y = -.2),
           autosize = TRUE
    )
}

# define UMAP output function for mod3 referring 'mt_plots_umap'
mod3_plots_umap <- function (D, title = "UMAP", 
                             ## scale argument
                             scale_data,  
                             ## color argument
                             color,
                             categorizing,
                             n_neighbors, 
                             ...) {
  X <- t(assay(D))
  if (any(is.na(X))) 
    stop("Data matrix for UMAP cannot contain NAs")
  ## dependent on scale checkbox
  if (scale_data) 
    X <- scale(X)
  um <- umap::umap(d = as.matrix(X), n_neighbors = as.numeric(n_neighbors))
  df <- data.frame(x = um$layout[, 1], y = um$layout[, 2], colData(D))
  colnames(df)[1:2] <- c("comp1", "comp2")
  ## reactivate plot title
  plot_title <- paste0(ifelse(scale_data,
                              "Scaled ", 
                              "Non-scaled "), title, 
                       " with ", n_neighbors, 
                       " neighbors colored by ", color)
  ## categorize coloring if color checkbox=T
  if(categorizing){
    df[, color] <- factor(df[, color])
  }
  # draw ggplot
  p <- ggplot(data = df,
              aes_string(x = "comp1", 
                         y = "comp2",
                         color=color)) + 
    geom_point() + 
    xlab("comp 1") + 
    ylab("comp 2") + 
    ggtitle(plot_title) +
    ## reactive legend title
    labs(color = color)
  
  # draw plotly
  ggplotly(p) %>% 
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = .5,
                         y = -.2),
           autosize = TRUE
    )
  
}
# load packages
library(shiny)
library(shinyWidgets)
library(maplet)
library(tidyverse)
library(DT)
library(plotly)

# Define UI for application
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
    # Six Head tabs to accommodate for navigation and comparison between modules
    tabPanel("Module 1", 
             sidebarLayout(
               sidebarPanel(id = "mod1_panel1",
                            # sidebar auto-scrolling with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>Module 1</b> requires extracting all the result objects one at a time."
                              )),
                            tags$p(
                              HTML("Users can assess results in a drop-down menu that offers a list of a stat_name and a plot type (e.g. “missingness”, “pval”)."
                              )),
                            br(),   
                            # select plot type or stats table
                            radioButtons("mod1_radio", "Select output type:",
                                         choices = list("Plot" = "plots", 
                                                        "Table" = "stats"),
                                         selected = "plots"
                            ),
                            br(),   
                            # define one UI object to select stat_name
                            uiOutput("mod1_select_statname_ui"),
                            br(),   
                            # define one UI object to select output type
                            uiOutput("mod1_select_object_ui"),
                            br(),
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                              )),
                            br(),
                            # delay the output
                            actionButton("mod1_go", "Update")
               ), 
               mainPanel(id = "mod1_panel2", 
                         # scrollable panel
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         br(), 
                         br(), 
                         br(), 
                         # dynamic number of plots
                         uiOutput('mod1_output')
               )
             )
    ), 
    tabPanel("Module 2", 
    ),
    tabPanel("Module 3", 
             sidebarLayout(
               sidebarPanel(id = "mod3_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>Module 3</b> requires generating an interactive 2D projection of PCA/UMAP."
                              )),
                            tags$p(
                              HTML("It displays a drop-down menu of all colData columns for coloring."
                              )),
                            # select one plot type
                            radioButtons("mod3_select_plot", "Select one plot type:", 
                                         choices = list("PCA" = "pca", 
                                                        "UMAP" = "umap")
                            ),
                            br(),   
                            # function argument
                            uiOutput("mod3_plot_argument"),
                            br(),
                            # select coloring colData and factor it
                            uiOutput("mod3_color_ui"),
                            br(),
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                              )),
                            br(),
                            # delay the output
                            actionButton("mod3_go", "Update")
               ), 
               mainPanel(id = "mod3_panel2", 
                         br(), 
                         br(), 
                         br(), 
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         # plotly
                         plotlyOutput('mod3_plot', height = 730)
               )
             )
    ),
    tabPanel("Module 4", "contents"),
    tabPanel("Module 5", "contents"),
    tabPanel("Module 6", "contents")
  )
)



# Define server logic required to draw outputs
server <- function(input, output) {
  ## create intermediate var to indicate coloring widgets
  inter_var <- reactive({
    if (input$mod3_select_plot=="pca" & input$mod3_pca_data_type=="scores") {
      "coloring"
    } else {
      "no_coloring"
    }
  })
  # create reactive inputs list
  mod1_input_object <- eventReactive(input$mod1_go, 
                                     {c(input$mod1_radio,
                                        input$mod1_select_statname,
                                        input$mod1_select_object)}
  )
  
  mod3_input_object <- eventReactive(input$mod3_go, 
                                     {c(input$mod3_select_plot, 
                                        input$mod3_select_colData,
                                        input$mod3_scale_data,
                                        input$mod3_checkbox_factor,
                                        input$mod3_pca_data_type,
                                        input$mod3_umap_n_neighbors)}
  )
  
  # create stat_name list dependent on radio button
  output$mod1_select_statname_ui <- renderUI({
    selectInput("mod1_select_statname", "Select one stat name:",
                width = "220px",
                choices = distinct(obj_name[obj_name$V1==input$mod1_radio, ], stat_name)$stat_name
    )
  })
  # create object list dependent on radio button and stat_name
  output$mod1_select_object_ui <- renderUI({
    selectInput("mod1_select_object", "Select one object:",
                width = "220px",
                choices = distinct(obj_name[obj_name$stat_name==input$mod1_select_statname & obj_name$V1==input$mod1_radio, ], V2)$V2
    )
  })
  # create reactive plotting argument for PCA/UMAP
  output$mod3_plot_argument <- renderUI({
    switch(
      input$mod3_select_plot,
      "pca"=list(selectInput("mod3_pca_data_type", "Select data type for PCA:",
                             width = "220px",
                             choices = c("scores", "loadings"),
                             selected = "scores"
      ),
      checkboxInput("mod3_scale_data", "Scaled data", 
                    value = TRUE
      )
      ),
      "umap"=list(numericInput("mod3_umap_n_neighbors", 
                               "Number of neighbors for UMAP:", 
                               value = 15,
                               width = "220px"
      ),
      checkboxInput("mod3_scale_data", "Scaled data", 
                    value = TRUE
      ),
      selectInput("mod3_select_colData", 
                  "Select one colData column:", 
                  choices = names(colData(D)),
                  selected = "BOX.NUMBER",
                  width = "220px"
      ),
      checkboxInput("mod3_checkbox_factor", 
                    "Categorical Coloring", 
                    value = FALSE
      )
      )
    )
  })
  output$mod3_color_ui <- renderUI({
    switch(
      inter_var(),
      "coloring"=list(selectInput("mod3_select_colData", 
                                  "Select one colData column:", 
                                  choices = names(colData(D)),
                                  selected = "BOX.NUMBER",
                                  width = "220px"
      ),
      checkboxInput("mod3_checkbox_factor", 
                    "Categorical Coloring", 
                    value = FALSE
      )
      ),
      "no_coloring"= NULL
    )
  })
  # Insert the right number of plot output objects into UI
  output$mod1_output_plot <- renderUI({
    ## limit plots to specified stat_name
    obj_name <- subset(obj_name, V1==mod1_input_object()[1])
    obj_name <- subset(obj_name, V2==mod1_input_object()[3])
    output_order <- obj_name %>%
      mutate(order=seq(from=1, to=n()))
    output_order <- subset(output_order, stat_name==mod1_input_object()[2])
    plots <- list()
    for(plot_i in seq_along(output_order$order)){
      plots[[plot_i]] <- mtm_res_get_entries(D, c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order[plot_i]]]
    }
    # there are multiple plots
    len_i <- length(plots)
    # some plots have multiple objects
    len_j <- length(plots[[1]]$output)
    # name every plot object in UI
    mod1_plot_output_list <- lapply(1:(len_i*len_j), function(i) {
      plotname <- paste("Plot", i, sep="")
      # locate the row in the `plots`
      row_n <- ceiling(i/len_j)
      ## set dynamic height of box scatter plots based on output2
      height <- if(plots[[1]]$fun[2]=="box"&plots[[1]]$fun[3]=="scatter"&!is.null(plots[[row_n]]$output2)){
        # (plots[[row_n]]$output2)*100
        1200
      } else if(plots[[1]]$fun[2]=="stats"){
        ## set dynamic height of plot_stats based on `nr` in output2
        (plots[[row_n]]$output2$nr)*20
      } else {
        560
      }
      plotOutput(plotname, height = height, width = 900)
    })
    # Convert the list to a tagList - this is necessary for the list of items
    # to display properly.
    do.call(tagList, mod1_plot_output_list)
  })
  
  # Call renderPlot for each one. Plots are only actually generated when they
  # are visible on the web page.
  for (i in 1:5) {
    # Need local so that each item gets its own number. Without it, the value
    # of i in the renderPlot() will be the same across all instances, because
    # of when the expression is evaluated.
    local({
      my_i <- i
      plotname <- paste("Plot", my_i, sep="")
      output[[plotname]] <- renderPlot({
        ## limit plots to specified stat_name
        obj_name <- subset(obj_name, V1==mod1_input_object()[1])
        obj_name <- subset(obj_name, V2==mod1_input_object()[3])
        output_order <- obj_name %>%
          mutate(order=seq(from=1, to=n()))
        output_order <- subset(output_order, stat_name==mod1_input_object()[2])
        plots <- list()
        for(plot_i in seq_along(output_order$order)){
          if(mod1_input_object()[1]=="plots"&mod1_input_object()[3]=="box"){
            ## only render top 10 pair objects for box-scatter plots
            dat <- mtm_res_get_entries(D, c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order[plot_i]]]$output
            dat[[1]]$data <- dat[[1]]$data %>% head(10000)
            plots[[plot_i]] <- dat
          } else {
            plots[[plot_i]] <- mtm_res_get_entries(D, c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order[plot_i]]]$output
          }
        }
        # there are multiple plots
        len_i <- length(plots)
        # some plots have multiple objects
        len_j <- length(plots[[1]])
        # locate the row in the `plots`
        row_n <- ceiling(my_i/len_j)
        # locate the column in the `plots`
        col_n <- ifelse((my_i %% len_j)==0, len_j, (my_i %% len_j))
        # render the plot object in each loop
        plots[[row_n]][col_n]
      })
    })
  }
  # render stats table of Mod1
  output$mod1_output_table <- renderDataTable({
    ## limit table to specified stat_name
    obj_name <- subset(obj_name, V1==mod1_input_object()[1])
    obj_name <- subset(obj_name, V2==mod1_input_object()[3])
    output_order <- obj_name %>%
      mutate(order=seq(from=1, to=n()))
    output_order <- subset(output_order, stat_name==mod1_input_object()[2])
    table <- mtm_res_get_entries(D, c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order]]$output$table
    
    ## put interested columns ahead
    table <- if ('term' %in% names(table)) {
      table %>%
        select(var, term, statistic, p.value, p.adj, everything())
    } else {
      table %>%
        select(var, statistic, p.value, p.adj, everything())
    }
    datatable(table,
              options = list(
                paging =TRUE,
                # limit number of rows
                pageLength =  10)) %>%
      # set two significant decimals
      formatSignif(columns=c("statistic", "p.value", "p.adj"), 
                  digits=2)
  })
  # render plots or table
  output$mod1_output <- renderUI({
    switch(
      mod1_input_object()[1],
      "plots" = uiOutput("mod1_output_plot"),
      "stats" = dataTableOutput("mod1_output_table")
    )
  })
  
  # render pca/umap of mod3
  output$mod3_plot <- renderPlotly({
    if (mod3_input_object()[1]=="pca"){
      mod3_plots_pca(D = D,
                     scale_data = mod3_input_object()[3],
                     color = mod3_input_object()[2],
                     categorizing=mod3_input_object()[4],
                     data_type = mod3_input_object()[5]
      )
    } else {
      mod3_plots_umap(D = D,
                      scale_data = mod3_input_object()[3],  
                      color = mod3_input_object()[2],
                      categorizing=mod3_input_object()[4],
                      n_neighbors = as.numeric(mod3_input_object()[6])
      )
    } 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)