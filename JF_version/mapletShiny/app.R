########################################################################## 
######################## run the result SE object ########################
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
######################## build the Shiny App #############################
########################################################################## 

D1 <- D

# Extract all the object names for accessor functions
obj_list <- data.frame()
for (i in seq_along(metadata(D1)$results)) {
  for (j in seq_along(metadata(D1)$results[[i]]$fun)) {
    obj_list[i, j] <- metadata(D1)$results[[i]]$fun[j]
  }
}

# Extract all the stat_name
stat_name <- data.frame(stat_name=NA)
for (i in seq_along(metadata(D1)$results)) {
  stat_name[i, 1] <- ifelse(is.null(metadata(D1)$results[[i]]$args$heading),
                            NA, metadata(D1)$results[[i]]$args$heading)
}
# distinct the not-null values
# stat_name <- distinct(subset(stat_name, !is.na(name)), name)

# merge object names and stat_name
obj_name <- cbind(obj_list, stat_name) %>%
  # impute stat_name with the most recent non-null value
  mutate(stat_name_fill=zoo::na.locf(stat_name, na.rm=FALSE))

## obj_name$name <- ifelse(is.na(obj_name$name), "(no stat_name)", obj_name$name)

# define plots & table objects
obj_type <- data.frame(plots=c("plots", NA), stats=c("stats", "post"))

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
                              condition = "input.mod1_radio!='plots'",
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
             absolutePanel(id = "mod3_panel1",
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
                           checkboxInput("mod3_checkbox_scale", "Scaled PCA", value = TRUE),
                           # check box for factor
                           checkboxInput("mod3_checkbox_factor", "Categorical Coloring", value = FALSE),
                           br(),
                           # dropdown menu for different assay
                           selectInput("mod3_select_assay", "Select one assay:", 
                                       choices = "Assays[1]",
                                       width = "220px"
                           ),
                           br(),   
                           # select one plot type
                           radioButtons("mod3_select_plot", "Select one plot type:", 
                                        choices = list("PCA" = "pca", 
                                                       "UMAP" = "umap")
                           ),
                           br(),   
                           # select one colData column for coloring
                           selectInput("mod3_select_colData", "Select one colData column:", 
                                       choices = names(colData(D)),
                                       selected = "BOX.NUMBER",
                                       width = "220px"
                           ),
                           br(),
                           tags$p(
                             HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                             )),
                           br(),
                           # delay the output
                           actionButton("mod3_go", "Update")
             ), 
             absolutePanel(id = "mod3_panel2", 
                           br(), 
                           br(), 
                           br(), 
                           style = "overflow-y: auto; position: absolute; left: 25%",
                           # plotly
                           plotlyOutput('mod3_plot', width = 1000, height = 700)
             )
    ),
    tabPanel("Module 4", "contents"),
    tabPanel("Module 5", "contents"),
    tabPanel("Module 6", "contents")
  )
)

# Define server logic required to draw outputs
server <- function(input, output) {
  # create stat_name list dependent on radio button
  output$mod1_select_statname_ui <- renderUI({
    selectInput("mod1_select_statname", "Select one stat name:",
                width = "220px",
                choices = distinct(obj_name[obj_name$V1==input$mod1_radio, ], stat_name_fill)$stat_name_fill
    )
  })
  # create object list dependent on radio button and stat_name
  output$mod1_select_object_ui <- renderUI({
    # create reactive object list
    mod1_radio_reactive <- reactive({obj_type[, input$mod1_radio]})
    selectInput("mod1_select_object", "Select one object:",
                width = "220px",
                choices = distinct(obj_name[obj_name$stat_name_fill==input$mod1_select_statname & obj_name$V1 %in% mod1_radio_reactive(), ], V2)$V2
    )
  })
  # create reactive output for plot
  mod1_output_object <- eventReactive(input$mod1_go, 
                                      {c(ifelse(input$mod1_select_object=="multtest", "post", input$mod1_radio), 
                                         input$mod1_select_object,
                                         input$mod1_select_statname)}
  )
  
  mod3_output_object <- eventReactive(input$mod3_go, 
                                      {c(input$mod3_select_plot, 
                                         input$mod3_select_colData,
                                         input$mod3_checkbox_scale,
                                         input$mod3_checkbox_factor)}
  )
  
  # define PCA output function for mod3 referring 'mt_plots_pca'
  mod3_plots_pca <- function(D, title = "PCA", 
                             ## scale argument
                             scale_data = mod3_output_object()[3],
                             ## color argument
                             color = mod3_output_object()[2],
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
    plot_title <- paste0(ifelse(mod3_output_object()[3],
                                "Scaled ", "Non-scaled "),
                         title, " - ",
                         mod3_output_object()[2])
    ## categorize coloring if color checkbox=T
    if(mod3_output_object()[4]){
      df[, mod3_output_object()[2]] <- factor(df[, mod3_output_object()[2]])
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
      labs(color = mod3_output_object()[2])
    
    # draw plotly
    ggplotly(p) %>% 
      layout(legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = .5,
                           y = -.2),
             autosize = TRUE
      )
  }
  
  # define UMAP output function for mod3 referring 'mt_plots_umap'
  mod3_plots_umap <- function (D, title = "UMAP", 
                               ## scale argument
                               scale_data = mod3_output_object()[3],  
                               ## color argument
                               color = mod3_output_object()[2],
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
    plot_title <- paste0(ifelse(mod3_output_object()[3],
                                "Scaled ", "Non-scaled "),
                         title, " - ",
                         mod3_output_object()[2])
    ## categorize coloring if color checkbox=T
    if(mod3_output_object()[4]){
      df[, mod3_output_object()[2]] <- factor(df[, mod3_output_object()[2]])
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
      labs(color = mod3_output_object()[2])
    
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
    ## limit the outputs to specified stat_name
    # output_order <- obj_name %>%
    #   filter(V1==mod1_output_object()[1], V2==mod1_output_object()[2]) %>%
    #   mutate(order=seq(from=1, to=n())) %>%
    #   filter(stat_name_fill==mod1_output_object()[3])
    obj_name <- subset(obj_name, V1==mod1_output_object()[1])
    obj_name <- subset(obj_name, V2==mod1_output_object()[2])
    output_order <- obj_name %>%
      mutate(order=seq(from=1, to=n()))
    output_order <- subset(output_order, stat_name_fill==mod1_output_object()[3])
    plots <- list()
    for(plot_i in seq_along(output_order$order)){
      plots[[plot_i]] <- mtm_res_get_entries(D1, c(mod1_output_object()[1], mod1_output_object()[2]))[[output_order$order[plot_i]]]$output
    }
    # there are multiple plots
    len_i <- length(plots)
    # some plots have multiple objects
    len_j <- length(plots[[1]])
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
        ## limit the outputs to specified stat_name
        # output_order <- obj_name %>%
        #   filter(V1==mod1_output_object()[1], V2==mod1_output_object()[2]) %>%
        #   mutate(order=seq(from=1, to=n())) %>%
        #   filter(stat_name_fill==mod1_output_object()[3])
        obj_name <- subset(obj_name, V1==mod1_output_object()[1])
        obj_name <- subset(obj_name, V2==mod1_output_object()[2])
        output_order <- obj_name %>%
          mutate(order=seq(from=1, to=n()))
        output_order <- subset(output_order, stat_name_fill==mod1_output_object()[3])
        plots <- list()
        for(plot_i in seq_along(output_order$order)){
          plots[[plot_i]] <- mtm_res_get_entries(D1, c(mod1_output_object()[1], mod1_output_object()[2]))[[output_order$order[plot_i]]]$output
        }
        # there are multiple plots
        len_i <- length(plots)
        # some plots have multiple objects
        len_j <- length(plots[[1]])
        # locate the row in the `plots`
        row_n <- ceiling(my_i/len_j)
        # locate the column in the `plots`
        col_n <- ifelse((my_i %% len_j)==0, len_j, (my_i %% len_j))
        
        plots[[row_n]][col_n]
      })
    })
  }
  # render stats table of Mod1
  output$mod1_table <- renderDataTable({
    ## limit the outputs to specified stat_name
    # output_order <- obj_name %>%
    #   filter(V1==mod1_output_object()[1], V2==mod1_output_object()[2]) %>%
    #   mutate(order=seq(from=1, to=n())) %>%
    #   filter(stat_name_fill==mod1_output_object()[3])
    obj_name <- subset(obj_name, V1==mod1_output_object()[1])
    obj_name <- subset(obj_name, V2==mod1_output_object()[2])
    output_order <- obj_name %>%
      mutate(order=seq(from=1, to=n()))
    output_order <- subset(output_order, stat_name_fill==mod1_output_object()[3])
    table <- mtm_res_get_entries(D1, c(mod1_output_object()[1], mod1_output_object()[2]))[[output_order$order]]$output$table
    datatable(table,
              options = list(
                paging =TRUE,
                # limit number of rows
                pageLength =  15)
    )
  })
  # render pca/umap of mod3
  output$mod3_plot <- renderPlotly({
    if (mod3_output_object()[1]=="pca"){
      mod3_plots_pca(D = D1)
    } else {
      mod3_plots_umap(D = D1)
    } 
  })
}

# Run the application 
shinyApp(ui = ui, server = server)