########################################################################## 
######################## run the result SE object ########################
########################################################################## 
library(maplet)
library(tidyverse)

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
######################## build the Shiny App #############################
########################################################################## 

# Extract all the object names for accessor functions
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
  }
}
# distinct the not-null values
# stat_name <- distinct(subset(stat_name, !is.na(name)), name)

# merge object names and stat_name
obj_name <- cbind(obj_list, stat_name)
obj_name$name <- ifelse(is.na(obj_name$name), "(no stat_name)", obj_name$name)


# load packages
library(shiny)
library(shinythemes)
library(plotly)
library(DT)

# Define UI for application
ui <- fluidPage(
  # remove shiny "red" warning messages on GUI
  tags$style(type="text/css",
             ".shiny-output-error { visibility: hidden; }",
             ".shiny-output-error:before { visibility: hidden; }"
  ),
  navbarPage(
    # embed Maplet logo and title
    title = div(img(src='logo.png',
                    style="margin-top: -10px; padding-right: 10px; padding-bottom: 10px",
                    height = 60),
                "BDS Capstone Maplet"),
    # sticky tabs while scrolling main panel
    position = c("fixed-top"),
    # set page theme
    theme = shinytheme("lumen"),
    # Six Head tabs to accommodate for navigation and comparison between modules
    tabPanel("Module 1", 
             absolutePanel(id = "mod1_panel1",
                           # sidebar auto-scrolling with main panel
                           style = "position:fixed; width: 20%;",
                           br(),   ## blank row
                           br(),   
                           br(),   
                           tags$p(
                             HTML("<b>Module 1</b> requires extracting all the result objects one at a time."
                             )),
                           tags$p(
                             HTML("Users can assess results in a drop-down menu that offers a list of a stat_name and a plot type (e.g. “missingness”, “pval”)."
                             )),
                           br(),   
                           # select one stat_name
                           selectInput("mod1_select_statname", "Select one stat name:", 
                                       choices = distinct(obj_name, name)$name,
                                       selected = "missingness",
                                       width = "220px"
                           ),
                           br(),   
                           # select plot type or stats table
                           radioButtons("mod1_radio", "Select output type:",
                                        choices = list("Plot" = "plots", 
                                                       "Table" = "stats"),
                                        selected = "plots"
                           ),
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
             conditionalPanel(id = "mod1_panel2", 
                              # Only show this panel if the plot type is selected
                              condition = "input.mod1_radio=='plots'",
                              # scrollable panel
                              style = "overflow-y: auto; position: absolute; left: 25%",
                              br(), 
                              br(), 
                              br(), 
                              # dynamic number of plots
                              uiOutput('mod1_plot')
             ),
             conditionalPanel(id = "mod1_panel3", 
                              # Only show this panel if the table type is selected
                              condition = "input.mod1_radio=='stats'",
                              # scrollable panel
                              style = "overflow-y: auto; position: absolute; left: 25%",
                              br(), 
                              br(), 
                              br(), 
                              dataTableOutput('mod1_table')
             )
    ), 
    tabPanel("Module 2", 
    ),
    tabPanel("Module 3", 
             absolutePanel(id = "mod2_panel1",
                           # sidebar autoscroll with main panel
                           style = "position:fixed; width: 20%;",
                           br(),   ## blank row
                           br(),   
                           br(),   
                           tags$p(
                             HTML("<b>Module 3</b> requires generating an interactive 2D projection of PCA/UMAP."
                             )),
                           tags$p(
                             HTML("It displays a drop-down menu of all colData columns for coloring."
                             )),
                           tags$p(
                             HTML("<b>Check boxes:</b>"
                             )),
                           # check box for scale
                           checkboxInput("mod2_checkbox_scale", "Scaled PCA", value = TRUE),
                           # check box for factor
                           checkboxInput("mod2_checkbox_factor", "Categorical Coloring", value = FALSE),
                           br(),
                           # dropdown menu for different assay
                           selectInput("mod2_select_assay", "Select one assay:", 
                                       choices = "",
                                       width = "220px"
                           ),
                           br(),   
                           # select one plot type
                           radioButtons("mod2_select_plot", "Select one plot type:", 
                                        choices = list("PCA" = "pca", 
                                                       "UMAP" = "umap")
                           ),
                           br(),   
                           # select one colData column for coloring
                           selectInput("mod2_select_colData", "Select one colData column:", 
                                       choices = names(colData(D1)),
                                       selected = "BOX.NUMBER",
                                       width = "220px"
                           ),
                           br(),
                           tags$p(
                             HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                             )),
                           br(),
                           # delay the output
                           actionButton("mod2_go", "Update")
             ), 
             absolutePanel(id = "mod2_panel2", 
                           br(), 
                           br(), 
                           br(), 
                           style = "overflow-y: auto; position: absolute; left: 25%",
                           # plotly
                           plotlyOutput('mod2_plot', width = 1000, height = 700)
             )
    ),
    tabPanel("Module 4", "contents"),
    tabPanel("Module 5", "contents"),
    tabPanel("Module 6", "contents")
  )
)

# Define server logic required to draw outputs
server <- function(input, output) {
  # create select input dependent on radio button and stat_name
  output$mod1_select_object_ui <- renderUI({
    selectInput("mod1_select_object", "Select one object:",
                width = "220px",
                choices = distinct(obj_name[obj_name$name==input$mod1_select_statname & obj_name$V1==input$mod1_radio, ], V2)$V2
    )
  })
  # create reactive output for plot
  mod1_output_object <- eventReactive(input$mod1_go, 
                                      {c(input$mod1_radio, input$mod1_select_object)}
  )
  
  mod2_output_object <- eventReactive(input$mod2_go, 
                                      {c(input$mod2_select_plot, 
                                         input$mod2_select_colData,
                                         input$mod2_checkbox_scale,
                                         input$mod2_checkbox_factor)}
  )
  
  # define PCA output function for Mod2 referring 'mt_plots_pca'
  mod2_plots_pca <- function(D, title = "PCA", 
                             ## scale argument
                             scale_data = mod2_output_object()[3],
                             ## color argument
                             color = mod2_output_object()[2],
                             pc1 = 1, pc2 = 2, 
                             data_type = "scores"){
    X = t(assay(D))
    if (any(is.na(X))) 
      stop("Data matrix for PCA cannot contain NAs")
    if (!(data_type %in% c("scores", "loadings"))) 
      stop("Show must be either 'scores' or 'loadings'")
    # scale scores if scale checkbox=T
    if (scale_data) 
      X <- scale(X)
    pca <- stats::prcomp(x = as.matrix(X), center = F, scale = F)
    expvar = (pca$sdev)^2/sum(pca$sdev^2)
    if (data_type == "scores") 
      df = data.frame(x = pca$x[, pc1], y = pca$x[, pc2], colData(D))
    else df = data.frame(x = pca$rotation[, pc1], y = pca$rotation[, 
                                                                   pc2], rowData(D))
    colnames(df)[1:2] <- c(sprintf("PC%d", pc1), sprintf("PC%d", 
                                                         pc2))
    ## reactivate axis labels
    pc1name <- sprintf("PC%d (%.1f%%)", pc1, expvar[pc1] * 
                         100)
    pc2name <- sprintf("PC%d (%.1f%%)", pc2, expvar[pc2] * 
                         100)
    ## reactivate plot title
    plot_title <- paste0(ifelse(mod2_output_object()[3],
                                "Scaled ", "Non-scaled "),
                         title, " - ",
                         mod2_output_object()[2])
    ## categorize coloring if color checkbox=T
    if(mod2_output_object()[4]){
      df[, mod2_output_object()[2]] <- factor(df[, mod2_output_object()[2]])
    }
    # draw ggplot
    p <- ggplot(data = df, 
                aes_string(
                  x = sprintf("PC%d", pc1), 
                  y = sprintf("PC%d", pc2),
                  color = color)
    ) + 
      geom_point() + 
      xlab(pc1name) + 
      ylab(pc2name) + 
      ggtitle(plot_title) +
      ## reactive legend title
      labs(color = mod2_output_object()[2])
    
    # draw plotly
    ggplotly(p) %>% 
      layout(legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = .5,
                           y = -.2),
             autosize = TRUE
      )
  }
  
  # define UMAP output function for Mod2 referring 'mt_plots_umap'
  mod2_plots_umap <- function (D, title = "UMAP", 
                               ## scale argument
                               scale_data = mod2_output_object()[3],  
                               ## color argument
                               color = mod2_output_object()[2],
                               n_neighbors = 15, 
                               ...) {
    X <- t(assay(D))
    if (any(is.na(X))) 
      stop("Data matrix for UMAP cannot contain NAs")
    ## dependent on scale checkbox
    if (scale_data) 
      X <- scale(X)
    um <- umap::umap(d = as.matrix(X), n_neighbors = n_neighbors)
    df <- data.frame(x = um$layout[, 1], y = um$layout[, 2], colData(D))
    colnames(df)[1:2] <- c("comp1", "comp2")
    ## reactivate plot title
    plot_title <- paste0(ifelse(mod2_output_object()[3],
                                "Scaled ", "Non-scaled "),
                         title, " - ",
                         mod2_output_object()[2])
    ## categorize coloring if color checkbox=T
    if(mod2_output_object()[4]){
      df[, mod2_output_object()[2]] <- factor(df[, mod2_output_object()[2]])
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
    labs(color = mod2_output_object()[2])
    
    # draw plotly
    ggplotly(p) %>% 
      layout(legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = .5,
                           y = -.2),
             autosize = TRUE
      )
    
  }
  
  # Insert the right number of plot output objects into UI
  output$mod1_plot <- renderUI({
    plots <- D1 %>% mtm_res_get_entries(mod1_output_object())
    # there are multiple plots
    len_i <- length(plots)
    # some plots have multiple objects
    len_j <- length(plots[[1]]$output)
    # name every plot object in UI
    mod1_plot_output_list <- lapply(1:(len_i*len_j), function(i) {
      plotname <- paste("Plot", i, sep="")
      plotOutput(plotname, height = 560, width = 900)
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
        plots <- D1 %>% mtm_res_get_entries(mod1_output_object())
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
  # render stats table of Mod1
  output$mod1_table <- renderDataTable({
    tables <- D1 %>% mtm_res_get_entries(mod1_output_object())
    datatable(tables[[1]]$output$table,
              options = list(
                paging =TRUE,
                # limit number of rows
                pageLength =  15)
    )
  })
  # render pca/umap of Mod2
  output$mod2_plot <- renderPlotly({
    if (mod2_output_object()[1]=="pca"){
      mod2_plots_pca(D = D1)
    } else {
      mod2_plots_umap(D = D1)
    } 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)
