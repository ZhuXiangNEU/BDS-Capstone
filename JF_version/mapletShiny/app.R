################################################################################
################################## Some Notes ##################################
################################################################################

# The codes are structured following some general notes of Shiny App, 
# you can just ignore this part if you are already familiar with Shiny. 

# Three steps to create a Shiny
# Step one: create layout of user interface (UI)
  # ui <- ...
# Step two: create rendering logic underlying the UI
  # server <- function(input, output) {}
# Step three: run the App
  # shinyApp(ui = ui, server = server)

# In ui, you generally define a sidebar panel for all the control widgets 
# and a main panel for the outputs.

# In the sidebar panel, commonly used widgets include
#   radioButtons() for single choice by clicking,
#   selectInput() for single choice by drop down box,
#   numericInput() for manually type in a number,
#   checkboxInput() for binary check box,
#   actionButton() for delayed outputs.
# You give each control widget an object_name which you refer to by input$object_name.

# In the main panel, you define output objects by 
#   plotOutput() for plot object,
#   plotlyOutput() for plotly object,
#   dataTableOutput() for table object,
#   downloadButton() for downloading output.
# You give each output object an object_name which you refer to by output$object_name.

# In server, you define rendering logic of UI objects by 
#   renderPlot() for plot objects,
#   renderPlotly() for plotly objects,
#   renderDataTable() for table objects.

# The input$object_name syntax can only be used inside of render(). 
# But you can refer to widget value outside of render() by creating reactive_object using reactive().

# In situations when some control widgets have dynamic values or 
# some outputs have dynamic structure, you define it as uiOutput().



################################################################################
######################### load packages and result SE object ###################
################################################################################

rm(list=ls())

# load packages
library(shiny)
library(shinyWidgets)
library(maplet)
# devtools::install_github(repo="krumsieklab/maplet@v1.0.1", subdir="maplet")
library(tidyverse)
library(DT)
library(plotly)
library(openxlsx)
library(readxl)

# load SE
load("SE.Rdata")


################################################################################
############ Create a data frame of object names and stat_name #################
##################### of all the result objects in SE ##########################
################################################################################

# Extract the object names
obj_list <- data.frame()
for (i in seq_along(metadata(D)$results)) {
  for (j in seq_along(metadata(D)$results[[i]]$fun)) {
    obj_list[i, j] <- metadata(D)$results[[i]]$fun[j]
  }
}

# Extract the stat_name
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

# Merge object names and stat_name
order_id <- 1:nrow(obj_list)
obj_name <- cbind(order_id, obj_list, stat_name) 
obj_name$stat_name <- ifelse(is.na(obj_name$stat_name), 
                             # name it as "(no stat_name)" for non-comparison analysis
                             "(no stat_name)", 
                             obj_name$stat_name)



################################################################################
######### Move some steps in advance to make server function codes #############
############## as few as possible to shorten waiting time ######################
######### as server function runs every time when control widget is updated ####
################################################################################

# assign an object of all stats tables
table_stats <- mtm_res_get_entries(D, c("stats", "univ"))
# assign an object of all stats plots
plot_stats <- mtm_res_get_entries(D, c("plots", "stats"))
# assign an object of all box plots
plot_box <- mtm_res_get_entries(D, c("plots", "box"))
# assign a data frame of all the object names of box plots
  # filter() cannot run in Shiny, use subset() instead
box_obj_name <- subset(obj_name, V1=="plots"&V2=="box")
box_output_order <- box_obj_name %>%
  mutate(order=seq(from=1, to=n()))

# define pathway annotation column (extracted from corresponding stat_bar
pwvar <- mtm_res_get_entries(D, c("plots", "stats"))[[1]]$args$group_col
# define threshold for significance (extracted from corresponding stat_bar plot)
alpha <- mtm_res_get_entries(D, c("plots", "stats"))[[1]]$args$feat_filter[[3]]

# get pathway annotations
rd <- rowData(D) %>% as.data.frame %>%
  dplyr::mutate(name=rownames(rowData(D))) %>%
  dplyr::select(name, !!sym(pwvar), BIOCHEMICAL)
if(class(rd[[pwvar]][1])=="AsIs"){
  # remove rows with missing pathway annotations and unnest pathway column
  miss_idx <- apply(rd, 1, function(x){x[[pwvar]][[1]] %>% is.null() %>% unname()})
  rd <- rd[!miss_idx,] %>% tidyr::unnest(!!sym(pwvar))
  # extract pathway_name column form pathways data frame
  rd %<>% dplyr::left_join(metadata(D)$pathways[[pwvar]], by=setNames("ID", pwvar)) %>% 
    dplyr::select(name, BIOCHEMICAL, pathway_name)
  # replace value for pathway column variable
  pwvar <- "pathway_name"
}

# reactive input with number of plots
n_plots <- 5
# function to remove a layer from a ggplot object
# (needed to strip the geom_text annotations)
# from https://stackoverflow.com/questions/13407236/remove-a-layer-from-a-ggplot2-chart
remove_geom <- function(ggplot2_object, geom_type) {
  # Delete layers that match the requested type.
  layers <- lapply(ggplot2_object$layers, function(x) {
    if (class(x$geom)[1] == geom_type) {
      NULL
    } else {
      x
    }
  })
  # Delete the unwanted layers.
  layers <- layers[!sapply(layers, is.null)]
  ggplot2_object$layers <- layers
  ggplot2_object
}
# access the plot data
x <- metadata(D)$results
x[[grep("plots_box_scatter", names(x), value = T)[2]]]$output[[1]]$data %<>%
  # add one column with order of plot to be displayed
  dplyr::mutate(rank=p.value %>% dplyr::dense_rank()) %>%
  # subselect only the first n_plots
  {.[.$rank<=n_plots,]} 
# plot reduced facets
x[[grep("plots_box_scatter", names(x), value = T)[2]]]$output[[1]] %>%
  remove_geom("GeomText")

## Get the data set 
get_data_by_name <- function(D, args.name, plot.nmae, stat.name) {
  plots <- mtm_res_get_entries(D, c("plots", plot.nmae))
  for (i in seq_along(plots)) {
    if (plots[[i]]$args[args.name] == stat.name) {
      data <- plots[[i]]$output[[1]]$data
    }
  }
  data
}

# create reverselog_trans for log10-SCALE in volcano plot
reverselog_trans <- function (base = exp(1)){
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf))
}

# define PCA output function for mod3 referring 'mt_plots_pca'
mod3_plots_pca <- function(D, title = "PCA", 
                           ## scale argument
                           scale_data,
                           ## color argument
                           color,
                           categorizing,
                           pc1 = 1, pc2 = 2, 
                           data_type,
                           hover,
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
    ## reactivate plot title
    plot_title <- paste0(ifelse(scale_data,
                                "Scaled ", 
                                "Non-scaled "),
                         title, " - ", color)
    ## categorize coloring if color checkbox=T
    if(categorizing){
      df[, color] <- factor(df[, color])
    }
    # customize hover text
    # hover_text <- paste0(names(data.frame(colData(D)))[as.numeric(hover)], ": ", 
    #                     data.frame(colData(D))[[as.numeric(hover)]])
    # draw ggplot
    p <- ggplot(data = df, 
                do.call(aes_string, as.list(structure(c("x","y",color,hover), names = c("x","y","colour",hover))))
                ) + 
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
    ## reactivate plot title
    plot_title <- paste0(ifelse(scale_data, 
                                "Scaled ", 
                                "Non-scaled "),
                         title)
    ## categorize coloring if color checkbox=T
    if(categorizing){
      df[, color] <- factor(df[, color])
    }
    # customize hover text
    # hover_text <- paste0(names(data.frame(rowData(D)))[as.numeric(hover)], ": ", 
    #                      data.frame(rowData(D))[[as.numeric(hover)]])
    # draw ggplot
    p <- ggplot(data = df, 
                do.call(aes_string, as.list(structure(c("x","y",color,hover), names = c("x","y","colour",hover))))
                ) + 
      geom_point() + 
      xlab(pc1name) + 
      ylab(pc2name) + 
      ggtitle(plot_title)
  }
  
  # draw plotly
  ggplotly(p, tooltip = c(color, hover)) %>%
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = 0.5, ## set position of legend
                         y = -0.2,
                         tracegroupgap = 5),
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
                             hover,
                             ...) {
  X <- t(assay(D))
  if (any(is.na(X))) 
    stop("Data matrix for UMAP cannot contain NAs")
  ## dependent on scale checkbox
  if (scale_data) 
    X <- scale(X)
  um <- umap::umap(d = as.matrix(X), n_neighbors = as.numeric(n_neighbors))
  df <- data.frame(x = um$layout[, 1], y = um$layout[, 2], colData(D))
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
  # customize hover text
  # hover_text <- paste0(names(data.frame(colData(D)))[as.numeric(hover)], ": ", 
  #                      data.frame(colData(D))[[as.numeric(hover)]])
  # draw ggplot
  p <- ggplot(data = df,
              do.call(aes_string, as.list(structure(c("x","y",color,hover), names = c("x","y","colour",hover))))
              ) + 
    geom_point() + 
    xlab("comp 1") + 
    ylab("comp 2") + 
    ggtitle(plot_title) +
    ## reactive legend title
    labs(color = color)
  
  # draw plotly
  ggplotly(p, tooltip = c(color, hover)) %>% 
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = .5,
                         y = -.2,
                         tracegroupgap = 5),
           autosize = TRUE
    )
  
}


################################################################################
########################## Define UI for Shiny application #####################
################################################################################

ui <- fluidPage(
  
  # set appearance customization -------------------------------------------------

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
    
  # Define layout of Module 1 ----------------------------------------------------

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
                                         selected = "stats"
                            ),
                            br(),   
                            # define one UI object to select stat_name
                            uiOutput("mod1_select_statname_ui"),
                            br(),   
                            # define one UI object to select output type
                            uiOutput("mod1_select_object_ui"),
                            br(),   
                            # how many bar plots to render
                            uiOutput("mod1_box_plot_num_ui"),
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

  # Define layout of Module 2 ----------------------------------------------------

  tabPanel("Module 2",
           sidebarLayout(
             sidebarPanel(
               id = "mod2_panel1",
               # sidebar auto-scrolling with main panel
               style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
               tags$p(
                 HTML(
                   "<b>Module 2:</b> StatsBar plot -> Equalizer/Volcano plot -> Box/Scatter plot."
                 )
               ),
               tags$p(
                 HTML(
                   "Given a SE and a statname, display a series of interactive plots at different granularities."
                 )
               ),
               br(),
               selectInput(
                 "mod2.stat",
                 "Select one stat name:",
                 choices = distinct(obj_name[obj_name$V1 == "plots" &
                                               obj_name$V2 == "stats",],
                                    stat_name)$stat_name
               ),
               br(),
               radioButtons(
                 "mod2.plot1",
                 "Select plot1 type:",
                 choices = list("Bar" = "bar",
                                "Not Bar" = "null"),
                 selected  = "bar"
               ),
               br(),
               radioButtons(
                 "mod2.plot2",
                 "Select plot2 type:",
                 choices = list("Equalizer" = "equalizer",
                                "Volcano" = "volcano"),
                 selected  = "volcano"
               ),
               br(),
               radioButtons(
                 "mod2.plot3",
                 "Select plot3 type:",
                 choices = list("Box" = "box",
                                "Scatter" = "scatter"),
                 selected  = "scatter"
               ),
               br(),
               checkboxInput("mod2.categorical", 
                             "Treat as categorical", 
                             value = FALSE
               ),
               br(),
               actionButton("mod2_go", "Update")
             ),
             mainPanel(
               id = "mod2_panel2",
               style = "overflow-y: auto; position: absolute; left: 25%",
               br(),
               br(),
               br(),
               # Bar plot or not
               uiOutput("mod2.p1"),
               br(),
               # equalizer or volcano
               uiOutput("mod2.p2"),
               br(),
               # box or scatter
               uiOutput("mod2.p3"),
               br()
             )
           )), 

  # Define layout of Module 3 ----------------------------------------------------

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
                            # function argument
                            uiOutput("mod3_pca_data"),
                            # select coloring colData and factor it
                            uiOutput("mod3_plot_argument"),
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
                         downloadButton("mod3_download_plotly", "download plotly"),
                         plotlyOutput('mod3_plot', height = 700)
               )
             )
    ),

  # Define layout of Module 4 --------------------------------------------------

  # Module 4
  tabPanel("Module 4",
           sidebarLayout(
             sidebarPanel(
               id = "mod4_panel1",
               # sidebar autoscroll with main panel
               style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
               tags$p(
                 HTML(
                   "<b>Module 4</b> requires collection on all statistical results in a table given one metabolite name."
                 )
               ),
               tags$p(
                 HTML(
                   "When clicking on one row, it should display interactive plots following the same orders in Module 2."
                 )
               ),
               # select one metabolite
               selectInput(
                 "mod4_metabolite",
                 "Select one metabolite:",
                 width = "220px",
                 choices = arrange(mtm_res_get_entries(D, c("stats", "univ"))[[1]]$output$table, var)$var,
                 selected = ""
               ),
               br(),
               checkboxInput("mod4.categorical", 
                             "Treat as categorical", 
                             value = FALSE
               ),
               br(),
               # delay the output
               actionButton("mod4_go", "Update")
             ),
             mainPanel(
               id = "mod4_panel2",
               br(),
               br(),
               br(),
               style = "overflow-y: auto; position: absolute; left: 25%",
               # stats table
               dataTableOutput('mod4_table'),
               br(),
               br(),
               # volcano plotly
               uiOutput("mod4.p1"),
               br(),
               br(),
               # box/scatter plotly
               uiOutput("mod4.p.ui"),
               uiOutput("mod4.p2")
             )
           )), 
  
# Define layout of Module 5 ----------------------------------------------------
    tabPanel("Module 5", 
             sidebarLayout(
               sidebarPanel(id = "mod5_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>Module 5</b> requires creating tables, distribution plots, or other graphics to explore the SE object."
                              )),
                            radioButtons("mod5_dimension", "Select one dimension:", 
                                         choices = list("Column Data" = "col", 
                                                        "Row Data" = "row")
                            ),
                            br(),
                            uiOutput("mod5_dimension_ui"),
                            br(),
                            # delay the output
                            actionButton("mod5_go", "Update")
               ), 
               mainPanel(id = "mod5_panel2", 
                         br(), 
                         br(), 
                         br(), 
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         downloadButton("mod5_download_plotly", "download plotly"),
                         plotlyOutput('mod5_plot', height = 600),
                         verbatimTextOutput("info"),
                         downloadButton("mod5_download_plotly2", "download plotly"),
                         plotlyOutput('mod5_plot2', height = 600)
               )
             )
    ),
# Define layout of Module 6 ----------------------------------------------------
    tabPanel("Module 6", 
             # Sidebar layout with input and output definitions ----
             sidebarLayout(
               # Sidebar panel for inputs ----
               sidebarPanel(
                 id = "mod6_panel1",
                 style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:absolute; width: 80%; overflow-y:auto; ",
                 # Input: Select a file ----
                 fileInput("file1", "Uploading File",
                           multiple = FALSE,
                           accept = c(".xlsx"),
                           width = "300px"),
                 # Input: Checkbox if file has header ----
                 checkboxInput("header", "Header", TRUE),
                 
                 tags$hr(),
                 
                 tags$p(HTML("<b>Sheets for Dimensions</b>")),
                 checkboxInput("mod6_assay_in_row", "Samples in rows?", TRUE),
                 tags$p(HTML("Assay sheet:")),
                 uiOutput("mod6_assay_sheet"),
                 tags$p(HTML("rowData sheet:")),
                 uiOutput("mod6_rowdata_sheet"),
                 tags$p(HTML("colData sheet:")),
                 uiOutput("mod6_coldata_sheet"),
                 actionButton("mod6_go", "Investigate"),
                 
                 tags$hr(),
                 
                 tags$p(HTML("<b>ID column for Dimensions</b>")),
                 tags$p(HTML("ID column in assay:")),
                 uiOutput("mod6_assay_id_column"),
                 tags$p(HTML("ID column in rowData:")),
                 uiOutput("mod6_rowdata_id_column"),
                 tags$p(HTML("ID column in colData:")),
                 uiOutput("mod6_coldata_id_column"),
                 
                 tags$hr(),
                 
                 tags$p(HTML("<b>Preprocessing</b>")),
                 tags$p(HTML("Max % missingness per feature:")),
                 numericInput("mod6_filter_feat_max", label = NULL,
                              value = 1,
                              min = 0,
                              max = 1,
                              step = 0.1,
                              width = "220px"),
                 tags$p(HTML("Max % missingness per feature (normalization):")),
                 numericInput("mod6_feat_max_norm", label = NULL,
                              value = 1,
                              min = 0,
                              max = 1,
                              step = 0.1,
                              width = "220px"),
                 tags$p(HTML("Max % missingness per sample:")),
                 numericInput("mod6_filter_sample_max", label = NULL,
                              value = 1,
                              min = 0,
                              max = 1,
                              step = 0.1,
                              width = "220px"),
                 tags$p(HTML("Sample coloring column:")),
                 uiOutput("mod6_pre_sample_color_column"),
                 tags$p(HTML("Batch column:")),
                 uiOutput("mod6_pre_batch_column"),
                 tags$p(HTML("PCA/UMAP coloring column:")),
                 uiOutput("mod6_pre_pca_color_column"),
                 tags$p(HTML("Heatmap annotation column:")),
                 uiOutput("mod6_pre_heatmap_anno_column"),
                 tags$p(HTML("Heatmap annotation row:")),
                 uiOutput("mod6_pre_heatmap_anno_row"),
                 
                 tags$hr(),
                 
                 tags$p(HTML("<b>Differential Analysis</b>")),
                 tags$p(HTML("Outcome variable:")),
                 uiOutput("mod6_outcome"),
                 checkboxInput("mod6_outcome_binary", "Binary outcome?", FALSE),
                 tags$p(HTML("Type of analysis:")),
                 selectInput("mod6_analysis_type", label = NULL,
                             width = "220px",
                             choices = c("lm","pearson","spearman","kendall"),
                             selected = "lm"),
                 tags$p(HTML("Multiple testing correction:")),
                 selectInput("mod6_mult_test_method", label = NULL,
                             width = "220px",
                             choices = c("BH","bonferroni","BY"),
                             selected = "BH"),
                 tags$p(HTML("Significance threshold:")),
                 numericInput("mod6_sig_threshold", label = NULL,
                              value = 0.05,
                              min = 0,
                              max = 1,
                              step = 0.01,
                              width = "220px"),
                 tags$p(HTML("Pathway aggregation in barplot:")),
                 uiOutput("mod6_group_col_barplot"),
                 tags$p(HTML("Barplot coloring column:")),
                 uiOutput("mod6_color_col_barplot"),
                 actionButton("mod6_go2", "LOAD RESULT SE")
               ),
                 
               # Main panel for displaying outputs ----
               mainPanel(
                 id = "mod6_panel2", 
                 style = "overflow-y: auto; position: absolute; left: 28%",
                 br(), 
                 br(), 
                 br(), 
                 # Output: Data file ----
                 downloadButton("download_se", "Download result SE .Rdata"),
                 br(),
                 br(),
                 dataTableOutput("mod6_assay"),
                 br(),
                 br(),
                 dataTableOutput("mod6_rowdata"),
                 br(),
                 br(),
                 dataTableOutput("mod6_coldata")
                 )
               )
             )
  )
)

################################################################################
################ Define server logic required to draw outputs ##################
################################################################################

server <- function(input, output) {

# Define rendering logic of control widgets in Module 1 ------------------------
  
  # create stat_name list dependent on radio button
  output$mod1_select_statname_ui <- renderUI({
    selectInput("mod1_select_statname", "Select one stat name:",
                width = "220px",
                choices = distinct(obj_name[obj_name$V1==input$mod1_radio, ], stat_name)$stat_name
    )
  })
  
  # create object list dependent on radio button and stat_name
  output$mod1_select_object_ui <- renderUI({
    if (input$mod1_radio=="stats"){
      NULL
    } else {
      selectInput("mod1_select_object", "Select one object:",
                width = "220px",
                choices = distinct(obj_name[obj_name$stat_name==input$mod1_select_statname&obj_name$V1==input$mod1_radio, ], V2)$V2
                )
    }
    
  })
  
  # create indicator of box plot output
  box_switch <- reactive({
    if (input$mod1_select_object=="box"){
      "box_plot"
    } else {
      "non_box_plot"
    }
  })
  
  ## get the order of selected stat_name
  ord <- reactive({
    if(input$mod1_select_statname %in% box_output_order$stat_name){
      box_output_order[box_output_order$stat_name==input$mod1_select_statname, ]$order
    } else {
      1
    }
  })
  
  ## create dynamic choices of number of box plots to render 
  output$mod1_box_plot_num_ui <- renderUI({
    switch(
      box_switch(),
      "box_plot"=list(selectInput("mod1_box_plot_num", 
                                  "Number of top box plots:", 
                                  choices = 1:(plot_box[[as.numeric(ord())]]$output2),
                                  width = "220px")
      ),
      "non_box_plot"=NULL
    )
  })
  
  # create reactive inputs list
  mod1_input_object <- eventReactive(input$mod1_go, ## delayed output
                                     {c(input$mod1_radio,
                                        input$mod1_select_statname,
                                        input$mod1_select_object,
                                        input$mod1_box_plot_num)}
                                     )
  
# Define rendering logic of outputs in Module 1 --------------------------------
  
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
        as.numeric(mod1_input_object()[4])*150
      } else {
        560
      }
      plotOutput(plotname, height = height, width = 850)
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
          plots[[plot_i]] <- mtm_res_get_entries(D, c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order[plot_i]]]$output
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
    table <- data.frame(var=row.names(rowData(D)), rowData(D)) %>%
      left_join(mtm_get_stat_by_name(D, mod1_input_object()[2]), 
                by=c("var"="var")
                ) %>%
      select(c(2, 20:26))
    
    ## put interested columns ahead
    table <- if ('term' %in% names(table)) {
      table %>%
        select(name, statistic, p.value, p.adj, term, everything()) %>%
        ## scientific notation
        mutate(statistic=formatC(statistic, format = "E", digits = 2),
               p.value=formatC(p.value, format = "E", digits = 2),
               p.adj=formatC(p.adj, format = "E", digits = 2),
               estimate=formatC(estimate, format = "E", digits = 2),
               std.error=formatC(std.error, format = "E", digits = 2)
        )
    } else {
      table %>%
        select(name, statistic, p.value, p.adj, everything())
    }
    datatable(table,
              options = list(
                # limit number of rows
                pageLength =  10,
                lengthMenu = c(10, 20, 50),
                ## set column width
                autoWidth = TRUE,
                columnDefs = list(list(width = '100px', targets = c(2:4))),
                scrollX = TRUE
              ))
  })
  # render plots or table
  output$mod1_output <- renderUI({
    switch(
      mod1_input_object()[1],
      "plots" = uiOutput("mod1_output_plot"),
      "stats" = dataTableOutput("mod1_output_table")
    )
  })

# Define rendering logic of outputs in Module 2 --------------------------------
  # Module 2: create reactive inputs list
  mod2_input_object <- eventReactive(input$mod2_go, 
                                     {c(input$mod2.stat,
                                        input$mod2.plot1,
                                        input$mod2.plot2,
                                        input$mod2.plot3,
                                        input$mod2.categorical)}
  )
  
  # Module 2: store reactive output plots
  session_store <- reactiveValues()
  
  # Module 2: plot 1
  output$mod2.p1 <- renderUI({
    inputs <- mod2_input_object()
    switch(
      inputs[2],
      "bar" = list(downloadButton("download_plotly_bar",
                                  "download bar plot"),
                   plotlyOutput("mod2.bar", height = 600)),
      "null" = NULL
    )
  })
  
  # Module 2: plot 1 - bar plot
  output$mod2.bar <- renderPlotly({
    inputs <- mod2_input_object()
    plots <- mtm_res_get_entries(D, c("plots", "stats"))
    for (i in seq_along(plots)) {
      if (plots[[i]]$args$stat_list == inputs[1]) {
        plot <- plots[[i]]$output[[1]]
      }
    }
    session_store$mod2.bar <- ggplotly(plot, source = "sub_bar") %>%
      layout(legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3))
    # render plotly graph
    session_store$mod2.bar
  })
  
  # Module 2: plot 1 - bar plot - html file
  output$download_plotly_bar <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # export plotly html widget as a temp file to download.
      saveWidget(as_widget(session_store$mod2.bar), file, selfcontained = TRUE)
    }
  )
  
  # Module 2: plot 2
  output$mod2.p2 <- renderUI({
    inputs <- mod2_input_object()
    d <- event_data("plotly_click", source = "sub_bar")
    
    vol_list <- list(
      downloadButton("download_plotly_volcano",
                     "download volcano plot"),
      plotlyOutput("mod2.vol", height = 600)
    )
    
    plot2 <- switch(inputs[3],
                    "equalizer" = switch(
                      inputs[2],
                      "bar" = if (!is.null(d)) {
                        list(
                          downloadButton("download_plotly_eq",
                                         "download equalizer plot"),
                          plotlyOutput("mod2.equal", height = 600)
                        )
                      },
                      "null" = list(
                        downloadButton("download_plotly_eq",
                                       "download equalizer plot"),
                        uiOutput("mod2.equal.ui"),
                        plotlyOutput("mod2.equal", height = 600)
                      )
                    ),
                    "volcano" = switch(inputs[2],
                                       "bar" = if (!is.null(d)) {
                                         vol_list
                                       },
                                       "null" = vol_list))
  })
  
  # Module 2: plot 2 - volcano plot
  output$mod2.vol <- renderPlotly({
    inputs <- mod2_input_object()
    
    # Get volcano data set if the plot1 is null
    data_vol <- get_data_by_name(D, "stat_name", "volcano", inputs[1])
    p.adj.significant <- alpha
    legend_name <- paste0("p.adj < ", p.adj.significant)
    
    # Get volcano data set if the plot1 is bar
    if (inputs[2] == "bar") {
      d <- event_data("plotly_click", source = "sub_bar")
      if (!is.null(d)) {
        # get the click information for the bar plot 
        data_bar <- get_data_by_name(D, "stat_list", "stats", inputs[1])
        lvls <- rev(levels(data_bar$label))
        label <- lvls[round(as.numeric(d$y))]
        sub_pathway_name <- data_bar[data_bar$label == label, ]$name
        # get the variable name by sub_pathway_name
        row_data <- rowData(D) %>% data.frame()
        names <- unlist(row_data[row_data[pwvar] == sub_pathway_name,]$name)
        data_vol <- data_vol[data_vol$name %in% names, ]
      }
    }
    # Set the legend color column
    data_vol[, legend_name] <- ifelse(data_vol$p.adj < as.numeric(p.adj.significant), TRUE, FALSE)
    
    # draw the plot2
    plot <- data_vol %>%
      ggplot(aes(x = statistic, y = p.value, color = !!sym(legend_name))) +
      geom_point() +
      scale_y_continuous(trans = reverselog_trans(10),
                         breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      labs(y = "p-value") +
      ggtitle(paste0(inputs[1]))
    
    if (inputs[2] == "bar") {
      plot <- plot +
        ggtitle(paste0(sub_pathway_name, "-", inputs[1]))
    }
    
    session_store$mod2.vol <- ggplotly(plot, source = "sub_vol") %>%
      layout(legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3))
    session_store$mod2.vol
  })
  
  # Module 2: plot 2 - volcano plot - html file
  output$download_plotly_volcano <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod2.vol), file, selfcontained = TRUE)
    }
  )
  
  # Module 2: plot 2 - equalizer plot - not bar
  output$mod2.equal.ui <- renderUI({
    inputs <- mod2_input_object()
    data_bar <- get_data_by_name(D, "stat_list", "stats", inputs[1])
    subpathways <- data_bar$name
    selectInput(
      "mod2.equal.path",
      "Select one pathway name:",
      choices = c(unique(unlist(subpathways))),
      selected = ""
    )
  })
  
  # Module 2: plot 2 - equalizer plot
  output$mod2.equal <- renderPlotly({
    inputs <- mod2_input_object()
    
    # add pathway annotations to results
    res <- maplet::mtm_get_stat_by_name(D, inputs[1]) %>%
      dplyr::left_join(rd, by=c("var"="name")) %>%
      dplyr::mutate(x = sign(statistic)*log10(p.adj)) %>%
      dplyr::filter(!is.na(BIOCHEMICAL))
    # compute multiple testing correction line
    sel <- res %>% 
      dplyr::filter(p.adj < alpha)
    if(nrow(sel)>0){
      xfine <- sel %>% .$p.adj %>% max(., na.rm = T)
    } else {
      xfine <- Inf
    }
    
    # If plot1 is not bar
    if (inputs[2] == "null") {
      df <- res %>%
        dplyr::filter(!!sym(pwvar)==input$mod2.equal.path)
    } else {
      d <- event_data("plotly_click", source = "sub_bar")
      if (!is.null(d)) {
        # get the click information for the bar plot
        data_bar <- get_data_by_name(D, "stat_list", "stats", inputs[1])
        lvls <- rev(levels(data_bar$label))
        label <- lvls[round(as.numeric(d$y))]
        sub_pathway_name <- data_bar[data_bar$label == label,]$name
        
        # create equalizer plots
        df <- res %>%
          dplyr::filter(!!dplyr::sym(pwvar) == sub_pathway_name)
      }
    }
    
    # colors
    clrs <- c("#9494FF","red")
    # x axis limits
    a = max(abs(df$x)) + 0.3
    plot <-
      ggplot(df, aes(x = x, y = BIOCHEMICAL)) +
      geom_vline(xintercept = 0, color = "gray") +
      (if (!is.infinite(xfine)) {
        geom_vline(
          xintercept = c(-log10(xfine), log10(xfine)),
          color = "red",
          alpha = 0.4
        )
      }) +
      (if (!is.infinite(xfine)) {
        ggtitle(sprintf("Differential Metabolites at alpha %.2f", alpha))
      } else{
        ggtitle(sprintf("No significant results at alpha %.2f", alpha))
      }) +
      geom_point(pch = 22,
                 fill = clrs[1],
                 size = 3) +
      facet_grid(as.formula(sprintf("%s~.", pwvar)),
                 scales = "free_y",
                 space = "free_y") +
      theme(
        strip.background = element_rect(fill = NA),
        strip.text = element_text(colour = 'black', face = "bold"),
        strip.text.y = element_text(angle = 0, hjust = 0),
        panel.grid.major.y = element_line(color = "gray"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_rect(fill = NA, color =
                                          "black")
      ) +
      ylab("") +
      xlab("sign(statistic)*log10(p.adj)") +
      scale_x_continuous(limits = c(-a, a))
    
    session_store$mod2.eq <- if (is.null(plot)) plotly_empty() else ggplotly(plot, source = "sub_eq")
    session_store$mod2.eq
  })
  
  # Module 2: plot 2 - equalizer plot - html file
  output$download_plotly_eq <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod2.eq), file, selfcontained = TRUE)
    }
  )
  
  # Module 2: plot 3 - box/scatter plot
  output$mod2.p3 <- renderUI({
    inputs <- mod2_input_object()
    d.eq <- event_data("plotly_click", source = "sub_eq")
    d.vol <- event_data("plotly_click", source = "sub_vol")
    
    download.name <- ifelse(inputs[4]=="box", "download box plot", "download scatter plot")
    plot.list <- list(
      downloadButton("download_plotly_box_scatter", download.name),
      plotOutput("mod2.box.scatter", height = 600)
    )
    
    if (!is.null(d.eq) | !is.null(d.vol))  {
      plot.list
    }
  })
  
  # Module 2: plot 3 - box/scatter plot
  output$mod2.box.scatter <- renderPlot({
    inputs <- mod2_input_object()
    # Get the data set
    data <- D %>%
      maplet:::mti_format_se_samplewise() %>%
      tidyr::gather(var, value, dplyr::one_of(rownames(D)))
    
    d.eq <- event_data("plotly_click", source = "sub_eq")
    d.vol <- event_data("plotly_click", source = "sub_vol")
    
    if (inputs[2] == "bar") {
      # get the click information for the bar plot
      d <- event_data("plotly_click", source = "sub_bar")
      
      if (!is.null(d)) {
        data_bar <- get_data_by_name(D, "stat_list", "stats", inputs[1])
        lvls <- rev(levels(data_bar$label))
        label <- lvls[round(as.numeric(d$y))]
        sub_pathway_name <- data_bar[data_bar$label == label,]$name
      }
    }
    
    # Get the metabolite name by click information
    if (inputs[3] == "equalizer") {
      if (!is.null(d.eq)) {
        res <- maplet::mtm_get_stat_by_name(D, inputs[1]) %>%
          dplyr::left_join(rd, by=c("var"="name")) %>%
          dplyr::mutate(x = sign(statistic)*log10(p.adj)) %>%
          dplyr::filter(!is.na(BIOCHEMICAL))
        
        # get the data set of the equalizer plot
        if (inputs[2] == "bar") {
          # Filter the data by selected pwvar
          df <- res %>%
            dplyr::filter(!!dplyr::sym(pwvar) == sub_pathway_name)
        } else {
          df <- res %>%
            dplyr::filter(!!dplyr::sym(pwvar) == input$mod2.equal.path)
        }
        # get the metabolite name by click info for equalizer plot
        metabolite <- df[as.numeric(d.eq$pointNumber) + 1,]$var[1]
        term <- df$term[1]
      }
    } else { # get the click info in volcano plot
      if (!is.null(d.vol)) {
        data_vol <- get_data_by_name(D, "stat_name", "volcano", inputs[1])
        if (inputs[2] == "bar") {
          # get the variable name by sub_pathway_name
          row_data <- rowData(D) %>% data.frame()
          names <- unlist(row_data[row_data[pwvar] == sub_pathway_name,]$name)
          data_vol <- data_vol[data_vol$name %in% names, ]
        }
        metabolite <- data_vol[as.numeric(d.vol$pointNumber) + 1,]$var[1]
        term <- data_vol$term[1]
        
      }
    }
    
    # Filter the data by metabolite name
    data <- data[data$var == metabolite, ]
    
    # Treat as categorical or not?
    if (inputs[5]) {
      data[, term] <- factor(data[, term])
    } else {
      data[, term] <- as.numeric(data[, term])
    }
    
    # Draw the plot3
    if (inputs[4] == "scatter") {
      plot <- data %>%
        ggplot(aes(x = !!sym(term), y = var)) +
        geom_point() +
        ggtitle(metabolite)
    } else {
      plot <- data %>%
        ggplot(aes(x = !!sym(term), y = var)) +
        geom_boxplot() +
        geom_jitter(width = 0.2) +
        ggtitle(metabolite)
    }
    
    session_store$mod2.box.scatter <- if (is.null(plot)) NULL else plot
    session_store$mod2.box.scatter
  })
  
  # Module 2: plot 3 - scatter/box plot - html file
  output$download_plotly_box_scatter <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file, plot = session_store$mod2.box.scatter, device = device)
    }
  )

# Define rendering logic of control widgets in Module 3 ------------------------
  output$mod3_pca_data <- renderUI({
    if(input$mod3_select_plot=="pca"){
      selectInput("mod3_pca_data_type", "Select data type for PCA:",
                  width = "220px",
                  choices = c("scores", "loadings"),
                  selected = "scores"
      )
    } else {
      NULL
    }
  })
  
  # create intermediate var to indicate coloring widgets
  inter_var <- reactive({
    if (input$mod3_select_plot=="pca" & input$mod3_pca_data_type=="scores") {
      "pca-scores"
    } else if(input$mod3_select_plot=="pca" & input$mod3_pca_data_type=="loadings"){
      "pca-loadings"
    } else {
      "umap"
    }
  })
  # create reactive plotting argument for PCA/UMAP
  output$mod3_plot_argument <- renderUI({
    switch(
      inter_var(),
      "pca-scores"=list(
      checkboxInput("mod3_scale_data", "Scaled data", 
                    value = TRUE
                    ),
      selectInput("mod3_select_colData", 
                  "Select one coloring variable:", 
                  choices = names(colData(D)),
                  selected = "BOX.NUMBER",
                  width = "220px"
                  ),
      checkboxInput("mod3_checkbox_factor", 
                    "Categorical Coloring", 
                    value = FALSE
                    ),
      selectInput("mod3_select_hover", 
                  "Select hovering text:", 
                  # selectInput coerces its output to character
                  # https://github.com/rstudio/shiny/issues/2367
                  # choices = setNames(seq_along(colData(D)), names(colData(D))),
                  choices = names(colData(D)),
                  selected = "sample",
                  width = "220px",
                  multiple=TRUE
                  )
      ),
      "pca-loadings"=list(
      checkboxInput("mod3_scale_data", "Scaled data", 
                    value = TRUE
                    ),
      selectInput("mod3_select_colData", 
                  "Select one coloring variable:", 
                  choices = names(rowData(D)),
                  selected = "SUPER_PATHWAY",
                  width = "220px"
                  ),
      checkboxInput("mod3_checkbox_factor", 
                    "Categorical Coloring", 
                    value = FALSE
                    ),
      selectInput("mod3_select_hover", 
                  "Select hovering text:", 
                  # choices = setNames(seq_along(rowData(D)), names(rowData(D))),
                  choices = names(rowData(D)),
                  selected = "name",
                  width = "220px",
                  multiple=TRUE
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
                  "Select one coloring variable:", 
                  choices = names(colData(D)),
                  selected = "BOX.NUMBER",
                  width = "220px"
                  ),
      checkboxInput("mod3_checkbox_factor", 
                    "Categorical Coloring", 
                    value = FALSE
                    ),
      selectInput("mod3_select_hover", 
                  "Select hovering text:", 
                  # choices = setNames(seq_along(colData(D)), names(colData(D))),
                  choices = names(colData(D)),
                  selected = "sample",
                  width = "220px",
                  multiple=TRUE
                  )
      )
    )
  })
  
  # create reactive inputs list
  mod3_input_object <- eventReactive(input$mod3_go, 
                                     {c(input$mod3_select_plot, 
                                        input$mod3_select_colData,
                                        input$mod3_scale_data,
                                        input$mod3_checkbox_factor,
                                        input$mod3_pca_data_type,
                                        input$mod3_umap_n_neighbors)}
  )

# Define rendering logic of outputs in Module 3 --------------------------------
  
  # render pca/umap of mod3
  output$mod3_plot <- renderPlotly({
    session_store$mod3_plotly <- if (mod3_input_object()[1]=="pca"){
      mod3_plots_pca(D = D,
                     scale_data = mod3_input_object()[3],
                     color = mod3_input_object()[2],
                     categorizing=mod3_input_object()[4],
                     data_type = mod3_input_object()[5],
                     hover = input$mod3_select_hover
      )
    } else {
      mod3_plots_umap(D = D,
                      scale_data = mod3_input_object()[3],  
                      color = mod3_input_object()[2],
                      categorizing=mod3_input_object()[4],
                      n_neighbors = as.numeric(mod3_input_object()[6]),
                      hover = input$mod3_select_hover
      )
    } 
    session_store$mod3_plotly
  })
  
  # download button
  output$mod3_download_plotly <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod3_plotly), file, selfcontained = TRUE)
    }
  )
  
# Define rendering logic of outputs in Module 4 --------------------------------
  
  # Module 4: general reactive stats table
  mod4_metabolite_table <-
    eventReactive(input$mod4_go,
                  {
                    table <- data.frame()
                    table_stats <- mtm_res_get_entries(D, c("stats", "univ"))
                    for (i in 2:length(table_stats)) {
                      tab <- table_stats[[i]]$output$table %>%
                        mutate(`stat name` = plot_stats[[i - 1]]$args$stat_list)
                      table <- rbind(table, tab)
                    }
                    table <- table %>%
                      select(var, statistic, p.value, p.adj, `stat name`, estimate, std.error) %>%
                      mutate(
                        statistic = formatC(statistic, format = "E", digits = 2),
                        p.value = formatC(p.value, format = "E", digits = 2),
                        p.adj = formatC(p.adj, format = "E", digits = 2),
                        estimate = formatC(estimate, format = "E", digits = 2),
                        std.error = formatC(std.error, format = "E", digits = 2)
                      ) %>%
                      filter(var == input$mod4_metabolite)
                  })
  
  # Module 4: output the stats table
  output$mod4_table <- renderDataTable({
    datatable(mod4_metabolite_table(),
              selection = "single",
              options = list(
                dom = 't',
                # limit number of rows
                pageLength =  10,
                lengthMenu = c(10, 20, 50)
              )
    )
  })  
  
  # mod4； catch the selected row
  observe({
    if (!is.null(input$mod4_table_rows_selected)) {
      session_store$mod4.tb.row <- input$mod4_table_rows_selected
    }
  })
  
  # mod4: extract the stat_name
  stat_name_selected <- reactive({
    mod4_metabolite_table() %>% 
      slice(round(as.numeric(session_store$mod4.tb.row))) %>%
      select(`stat name`)
  })
  
  # Module 4: volcano plot
  output$mod4.p1 <- renderUI({
    if (!is.null(session_store$mod4.tb.row)) {
      list(
        downloadButton("mod4_download_plotly_volcano", "download volcano plot"),
        plotlyOutput('mod4_volcano', height = 800)
      )
    }
  })
  
  # Module 4: volcano plot by using stat_name
  output$mod4_volcano <- renderPlotly({
    # Get volcano data set
    data_vol <- get_data_by_name(D, "stat_name", "volcano", stat_name_selected())
    isSelected <- input$mod4_metabolite
    
    # Set the legend color column
    data_vol[, "isSelected"] <- ifelse(data_vol$var==isSelected, TRUE, FALSE)
    
    plot <- data_vol %>%
      ggplot(aes(x = statistic, y = p.value, color = isSelected)) +
      geom_point() +
      scale_y_continuous(trans = reverselog_trans(10),
                         breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      labs(y = "p-value") +
      ggtitle(paste0(stat_name_selected()))
    
    session_store$mod4.vol <- ggplotly(plot, source = "mod4_sub_vol") %>%
      layout(legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3))
    session_store$mod4.vol
  })
  
  # Module 4: volcano plot - html file
  output$mod4_download_plotly_volcano <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod4.vol), file, selfcontained = TRUE)
    }
  )
  
  # Module 4: box/scatter plot
  output$mod4.p2 <- renderUI({
    d <- event_data("plotly_click", source = "mod4_sub_vol")
    
    if (!is.null(d)) {
      download.name <- ifelse(
        input$mod4.box.or.scatter == "box",
        "download box plot",
        "download scatter plot"
      )
      list(
        downloadButton("mod4_download_box_scatter", download.name),
        plotOutput("mod4.box.scatter", height = 600)
      )
    }
  })
  
  # Module 4: box/scatter - ui
  output$mod4.p.ui <- renderUI({
    d <- event_data("plotly_click", source = "mod4_sub_vol")
    if (!is.null(d)) {
      radioButtons(
        "mod4.box.or.scatter",
        "Select plot type:",
        choices = list("Box" = "box",
                       "Scatter" = "scatter"),
        selected  = "scatter"
      )
    }
  })
  
  # Module 4: box/scatter plot
  output$mod4.box.scatter <- renderPlot({
    # Get the data set
    data <- D %>%
      maplet:::mti_format_se_samplewise() %>%
      tidyr::gather(var, value, dplyr::one_of(rownames(D)))
    
    d <- event_data("plotly_click", source = "mod4_sub_vol")
    
    if (!is.null(d)) {
      data_vol <- get_data_by_name(D, "stat_name", "volcano", stat_name_selected())
      metabolite <- data_vol[as.numeric(d$pointNumber) + 1,]$var[1]
      term <- data_vol$term[1]
      
      # Filter the data by metabolite name
      data <- data[data$var == metabolite, ]
      
      # Treat as categorical or not?
      if (input$mod4.categorical) {
        data[, term] <- factor(data[, term])
      } else {
        data[, term] <- as.numeric(data[, term])
      }
      
      # Draw the plot
      if (input$mod4.box.or.scatter == "scatter") {
        plot <- data %>%
          ggplot(aes(x = !!sym(term), y = var)) +
          geom_point() +
          ggtitle(metabolite)
      } else {
        plot <- data %>%
          ggplot(aes(x = !!sym(term), y = var)) +
          geom_boxplot() +
          geom_jitter(width = 0.2) +
          ggtitle(metabolite)
      }
    }
    
    session_store$mod4.box.scatter <- if (is.null(plot)) NULL else plot
    session_store$mod4.box.scatter
  })
  
  # Module 4: scatter/box plot - png file
  output$mod4_download_box_scatter <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".png", sep = "")
    },
    content = function(file) {
      device <- function(..., width, height) grDevices::png(..., width = width, height = height, res = 300, units = "in")
      ggsave(file, plot = session_store$mod4.box.scatter, device = device)
    }
  )
  
# Define rendering logic of control widgets in Module 5 ----------------------
  output$mod5_dimension_ui <- renderUI({
    switch(input$mod5_dimension,
           "col"=list(selectInput("mod5_var1_select", 
                                  "Select the primary variable:", 
                                  choices = names(colData(D)),
                                  selected = "Age",
                                  width = "220px"),
                      checkboxInput("mod5_var1_type", 
                                    "Continuous", 
                                    value = TRUE),
                      selectInput("mod5_var2_select", 
                                  "Select the primary variable:", 
                                  choices = names(colData(D)),
                                  selected = "sample",
                                  width = "220px"),
                      checkboxInput("mod5_var2_type", 
                                    "Continuous", 
                                    value = TRUE)
                      ),
           "row"=selectInput("mod5_rowdata_plot", 
                             "Select one plot for row data:", 
                             choices = c("SUPER_PATHWAY"),
                             width = "220px")
    )
  })
  
# Define rendering logic of outputs in Module 5 ------------------------------
  mod5_input <- eventReactive(input$mod5_go,{
    c(input$mod5_var1_select,
      input$mod5_var1_type,
      input$mod5_var2_select,
      input$mod5_var2_type,
      input$mod5_rowdata_plot)
  })
  
  output$mod5_plot <- renderPlotly({
    session_store$mod5_plotly <- switch(input$mod5_dimension,
           "col"=
             if(mod5_input()[2]==TRUE & mod5_input()[4]==TRUE){
      p <- ggplot(as.data.frame(colData(D)),
             aes(!!sym(mod5_input()[3]), !!sym(mod5_input()[1]))
             ) + geom_point()
      ggplotly(p)
    } else if(mod5_input()[2]==TRUE & mod5_input()[4]==FALSE) {
      p <- ggplot(as.data.frame(colData(D)),
             aes(!!sym(mod5_input()[3]), !!sym(mod5_input()[1]))
      ) + geom_boxplot(outlier.shape = NA) + geom_jitter(width=.1)
      ggplotly(p)
    } else if(mod5_input()[2]==FALSE & mod5_input()[4]==TRUE) {
      p <- ggplot(as.data.frame(colData(D)),
             aes(!!sym(mod5_input()[1]), !!sym(mod5_input()[3]))
      ) + geom_boxplot(outlier.shape = NA) + geom_jitter(width=.1)
      ggplotly(p)
    } else {
      p <- ggplot(as.data.frame(colData(D)),
             aes(x=!!sym(mod5_input()[3]), fill=!!sym(mod5_input()[1]))
      ) + geom_bar()
      ggplotly(p)
    },
           "row"=
      as.data.frame(rowData(D)) %>%
      rename(var=mod5_input()[5]) %>%
      group_by(var) %>%
      summarise(count=n()) %>%
      plot_ly(labels = ~var, values = ~count, type = 'pie',
              source="mod5-click") %>% 
      layout(legend = list(orientation = "h",   # show entries horizontally
                           xanchor = "center",  # use center of legend as anchor
                           x = .5,
                           y = -.2,
                           tracegroupgap = 5),
             autosize = TRUE
      )
    )
    session_store$mod5_plotly
    }
  )
  # download button
  output$mod5_download_plotly <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod5_plotly), file, selfcontained = TRUE)
    }
  )
  
  output$info <- renderPrint({
    d5 <- event_data("plotly_click", source = "mod5-click")
    if(!is.null(d5)){
      d5
    }
  })
  
  output$mod5_plot2 <- renderPlotly({
    d5 <- event_data("plotly_click", source = "mod5-click")
    pie_dat <- as.data.frame(rowData(D))
    
    if (!is.null(d5)){
      lvls <- rev(pie_dat$SUPER_PATHWAY)
      label <- lvls[round(as.numeric(d5$pointNumber))+1]
      
      session_store$mod5_plot2 <- 
        pie_dat[pie_dat$SUPER_PATHWAY == label, ] %>%
        rename(var="SUB_PATHWAY") %>%
        group_by(var) %>%
        summarise(count=n()) %>%
        plot_ly(labels = ~var, values = ~count, type = 'pie') %>% 
        layout(legend = list(orientation = "h",   # show entries horizontally
                             xanchor = "center",  # use center of legend as anchor
                             x = .5,
                             y = -.2,
                             tracegroupgap = 5),
               autosize = TRUE)
      session_store$mod5_plot2
    }
  })
  # download button
  output$mod5_download_plotly2 <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod5_plotly2), file, selfcontained = TRUE)
    }
  )
  
  # Define rendering logic of control widgets in Module 6 ----------------------
  
  output$mod6_assay_sheet <- renderUI({
    req(input$file1)
    selectInput("assay_sheet", label = NULL,
                width = "220px",
                choices = getSheetNames(as.character(input$file1$datapath))
    )
  })
  df_assay <- reactive({
    read_excel(as.character(input$file1$datapath),
               col_names = input$header,
               sheet=input$assay_sheet)
    })
  
  output$mod6_rowdata_sheet <- renderUI({
    req(input$file1)
    selectInput("rowdata_sheet", label = NULL,
                width = "220px",
                choices = getSheetNames(as.character(input$file1$datapath))
                )
  })
  df_rowdata <- reactive({
    read_excel(as.character(input$file1$datapath),
               col_names = input$header,
               sheet=input$rowdata_sheet)
  })
  
  output$mod6_coldata_sheet <- renderUI({
    req(input$file1)
    selectInput("coldata_sheet", label = NULL,
                width = "220px",
                choices = getSheetNames(as.character(input$file1$datapath))
    )
  })
  df_coldata <- reactive({
    read_excel(as.character(input$file1$datapath),
               col_names = input$header,
               sheet=input$coldata_sheet)
  })
  
  output$mod6_assay_id_column <- renderUI({
    selectInput("assay_id_column", label = NULL,
                width = "220px",
                choices = colnames(df_assay())
    )
  })
  
  output$mod6_rowdata_id_column <- renderUI({
    selectInput("rowdata_id_column", label = NULL,
                width = "220px",
                choices = colnames(df_rowdata())
    )
  })
  
  output$mod6_coldata_id_column <- renderUI({
    selectInput("coldata_id_column", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_sample_color_column <- renderUI({
    selectInput("pre_sample_color_column", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_batch_column <- renderUI({
    selectInput("pre_batch_column", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
 
  output$mod6_pre_pca_color_column <- renderUI({
    selectInput("pre_pca_color_column", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_heatmap_anno_column <- renderUI({
    selectInput("pre_heatmap_anno_column", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_heatmap_anno_row <- renderUI({
    selectInput("pre_heatmap_anno_row", label = NULL,
                width = "220px",
                choices = colnames(df_rowdata())
    )
  })
  
  output$mod6_outcome <- renderUI({
    selectInput("outcome", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_group_col_barplot <- renderUI({
    selectInput("group_col_barplot", label = NULL,
                width = "220px",
                choices = colnames(df_rowdata())
    )
  })
  
  output$mod6_color_col_barplot <- renderUI({
    selectInput("color_col_barplot", label = NULL,
                width = "220px",
                choices = colnames(df_rowdata())
    )
  })
  
  # Define rendering logic of outputs in Module 6 ------------------------------
  mod6_filepath <- 
    eventReactive(input$mod6_go, ## delayed output
                  {c(input$file1$datapath)
                    })
  
  output$mod6_assay <- renderDataTable({
    table <- read_excel(as.character(mod6_filepath()),
                        col_names = input$header,
                        sheet=input$assay_sheet)
    datatable(table,
              caption="Original Assay Data",
              options = list(
                # limit number of rows
                pageLength =  10,
                lengthMenu = c(10, 20, 50),
                autoWidth = TRUE
              ))
  })
  
  output$mod6_rowdata <- renderDataTable({
    table <- read_excel(as.character(mod6_filepath()),
                        col_names = input$header,
                        sheet=input$rowdata_sheet)
    datatable(table,
              caption="Original rowData",
              options = list(
                # limit number of rows
                pageLength =  10,
                lengthMenu = c(10, 20, 50),
                autoWidth = TRUE
              ))
  })
  
  output$mod6_coldata <- renderDataTable({
    table <- read_excel(as.character(mod6_filepath()),
                        col_names = input$header,
                        sheet=input$coldata_sheet)
    datatable(table,
              caption="Original colData",
              options = list(
                # limit number of rows
                pageLength =  10,
                lengthMenu = c(10, 20, 50),
                autoWidth = TRUE
              ))
  })
  
  # https://mastering-shiny.org/action-transfer.html
  output$download_se <- downloadHandler(
    filename = function() {
      paste0("SE_", Sys.Date(), ".Rdata")
    },
    content = function(file) {
      ## Loading Data ----
      # file for assay
      file_data <- "Module6_Data.xlsx"
      # data sheet
      sheet_data <- "data" # to be choosen by dropdown menu displaying available sheet names from file_data
      # samples in rows?
      sir <- F # checkbox for T/F values
      # id column
      id_col_data <- "BIOCHEMICAL" # to be chosen by dropdown menu from column names
      # file for rowData
      file_feature <- file_data # should be the same as data_file by default, with the option to select a different file
      # data sheet
      sheet_feature <- "metanno" # to be choosen by dropdown menu from available sheet names from file_feature
      # id column
      id_col_feature <- "BIOCHEMICAL" # to be chosen by dropdown menu from column names
      # file for coData
      file_sample <- file_data # should be the same as data_file by default, with the option to select a different file
      # data sheet
      sheet_sample <- "sampleanno" # to be choosen by dropdown menu from available sheet names from file_feature
      # id column
      id_col_sample <- "SAMPLE_NAME" # to be chosen by dropdown menu from column names
      ## Preprocessing ----
      # maximal allowed missingness percentage per feature
      filter_feat_max <- 0.5 # can be any number between 0 and 1. Can be left empty. Default 1.
      # maximal allowed missingness percentage per sample
      filter_sample_max <- 1 # can be any number between 0 and 1. Can be left empty. Default 1.
      # color column for sample boxplot
      sample_boxplot_col <- "Tissue"
      # max missingness percentage per feature for normalization
      feat_max_norm <- 0.2 # can be any number between 0 and 1. Can be left empty. Default 1.
      # batch column
      batch <- "Batch_RUN_DAY" # to be chosen by dropdown menu from colData columns. Can be left empty. Default NULL.
      # columns for PCA/UMAP coloring
      pca_cols <- c("Batch_RUN_DAY","Tissue","Gleason") # to be chosen by dropdown menu from colData columns. Multiple choice possible. Default NULL.
      # columns for heatmap annotations
      heatmap_anno_col <- pca_cols # to be chosen by dropdown menu from colData columns. Multiple choice possible. Default NULL.
      heatmap_anno_row <- "SUPER_PATHWAY" # to be chosen by dropdown menu from rowData columns. Multiple choice possible. Default NULL.
      ## Differential Analysis ----
      # outcome
      var="Gleason" # to be chosen by dropdown menu from colData columns.
      # is the outcome binary?
      binary=F 
      # type of analysis
      analysis_type="kendall" # to be chosen from dropdown menu between the following: c("lm","pearson","spearman","kendall"). Default "lm".
      # method for multiple testing correction
      mult_test_method="BH" # to be chosen from dropdown menu between the following: c("BH","bonferroni","BY"). Default "BH".
      # significance threshold
      alpha=0.05 # can be any number between 0 and 1. Default 0.05.
      # column for pathway aggregation in barplot
      group_col_barplot="SUB_PATHWAY" # to be chosen by dropdown menu from rowData columns.
      # column for barplot colors
      color_col_barplot="SUPER_PATHWAY" # to be chosen by dropdown menu from rowData columns. Can be left empty. Default NULL.
      # Define Functions ----
      diff_analysis_func <- function(D,
                                     var,
                                     binary=F,
                                     analysis_type="lm",
                                     mult_test_method="BH",
                                     alpha=0.05,
                                     group_col_barplot,
                                     color_col_barplot=NULL){
        D %<>%
          mt_reporting_heading(heading = sprintf("%s Differential Analysis",var), lvl = 2) %>%
          {.}
        
        if(analysis_type=="lm"){
          D %<>%
            mt_stats_univ_lm(formula = as.formula(sprintf("~  %s",var)), stat_name = sprintf("%s analysis",var)) %>%
            {.}
        }else{
          D %<>%
            mt_stats_univ_cor(in_col = var, stat_name = sprintf("%s analysis",var),method = analysis_type) %>%
            {.}
        }
        
        if(binary){
          D %<>%
            mt_post_fold_change(stat_name = sprintf("%s analysis",var))
        }
        D %<>%
          mt_post_multtest(stat_name = sprintf("%s analysis",var), method = mult_test_method) %>%
          mt_reporting_stats(stat_name = sprintf("%s analysis",var), stat_filter = p.adj < alpha) %>%
          mt_plots_volcano(stat_name = sprintf("%s analysis",var),
                           x = !!sym(ifelse(binary,"fc","statistic")),
                           feat_filter = p.adj < alpha,
                           color = p.adj < alpha) %>%
          mt_plots_box_scatter(stat_name = sprintf("%s analysis",var),
                               x = !!sym(var),
                               plot_type = ifelse(binary,"box","scatter"),
                               feat_filter = p.adj < alpha, 
                               feat_sort = p.value,
                               annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}") %>%
          mt_plots_stats_pathway_bar(stat_list = sprintf("%s analysis",var),
                                     y_scale = "count",
                                     feat_filter = p.adj < alpha,
                                     group_col = group_col_barplot,
                                     color_col = color_col_barplot) %>%
          {.}
        
        return(D)
      }
      # Loading Data ----
      D <-
        mt_load_xls(file=file_data, sheet=sheet_data, samples_in_row=sir, id_col=id_col_data) %>%
        mt_anno_xls(file=file_feature, sheet=sheet_feature, anno_type="features", anno_id_col=id_col_feature, data_id_col = "name") %>%
        mt_anno_xls(file=file_sample, sheet=sheet_sample, anno_type="samples", anno_id_col =id_col_sample, data_id_col ="sample") %>%
        mt_reporting_data() %>%
        {.}
      # Preprocessing ----
      D <- D %>%  
        mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
        
        mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
        mt_plots_missingness(feat_max=filter_feat_max,samp_max = filter_sample_max) %>%
        mt_pre_filter_missingness(feat_max = filter_feat_max, samp_max = filter_sample_max) %>%
        mt_plots_missingness(feat_max=filter_feat_max, samp_max = filter_sample_max) %>%
        mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
        mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
        
        mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
        mt_plots_sample_boxplot(color=!!sym(sample_boxplot_col), title = "Original", plot_logged = T) %>%
        {.}
      if(!is.null(batch)){
        D %<>%
          mt_pre_batch_median(batch_col = batch)
      }
      D <- D %>%
        mt_plots_sample_boxplot(color=!!sym(sample_boxplot_col), title = "After batch correction", plot_logged = T) %>%
        mt_pre_norm_quot(feat_max = feat_max_norm) %>%
        mt_plots_dilution_factor(in_col=sample_boxplot_col) %>%
        mt_plots_sample_boxplot(color=!!sym(sample_boxplot_col), title = "After normalization", plot_logged = T) %>%
        mt_pre_trans_log() %>%
        mt_pre_impute_knn() %>%
        mt_plots_sample_boxplot(color=!!sym(sample_boxplot_col), title = "After imputation", plot_logged = T) %>%
        mt_pre_outlier_detection_univariate() %>%
        mt_reporting_data() %>%
        
        # mt_reporting_heading(heading = "Get Pathway Annotations", lvl = 1) %>%
        # mt_anno_hmdb_to_kegg(in_col = "HMDb", out_col = "KEGG_ids") %>%
        # mt_anno_pathways_hmdb(in_col = "HMDb", out_col = "pathway", pwdb_name = "KEGG") %>%
        # mt_anno_pathways_remove_redundant(feat_col = "KEGG_ids", pw_col = "pathway") %>%
        
        mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
        {.}
      # add PCA/UMAP plots
      lapply(pca_cols, function(x){
        D <<- D %>%
          mt_plots_pca(scale_data = T, title = sprintf("scaled PCA - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
          mt_plots_umap(scale_data = T, title = sprintf("scaled UMAP - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
          {.}
      }) %>% invisible
      # add heatmap
      D %<>%
        mt_plots_heatmap(scale_data = T, annotation_col = heatmap_anno_col, annotation_row = heatmap_anno_row,
                         clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
        {.}
      # Differential analysis ----
      D %<>%
        mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
        
        diff_analysis_func(var=var,binary=binary,analysis_type=analysis_type, mult_test_method=mult_test_method,alpha=alpha,
                           group_col_barplot=group_col_barplot,color_col_barplot=color_col_barplot) %>%
        {.}
      # write to local
      saveRDS(D, file)
    }
  )
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)