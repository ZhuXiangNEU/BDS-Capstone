rm(list=ls())

################################################################################
####################load packages, functions and result SE object ##############
################################################################################

# load packages
# devtools::install_github(repo="krumsieklab/maplet@v1.0.1", subdir="maplet")
library(shiny)
library(shinyWidgets)
library(maplet)
library(shinydashboard)
library(shinyjs)
library(tidyverse)
library(DT)
library(plotly)
library(openxlsx)
library(readxl)
library(RColorBrewer)
# refer help functions
source("help_functions.R")

# load SE with fixed name
load("SE.Rdata")
# extract object names from result SE 'D'
obj_name <- get_obj_name(D)
# define pathway annotation column (extracted from corresponding stat_bar
pwvar <- mtm_res_get_entries(D, c("plots", "stats"))[[1]]$args$group_col
# define threshold for significance (extracted from corresponding stat_bar plot)
alpha <- mtm_res_get_entries(D, c("plots", "stats"))[[1]]$args$feat_filter[[3]]
# get pathway annotations
rd <- get_pathway_annotations(D, pwvar)


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
  # adjust tab height
  tags$head(
    tags$style(HTML(' .navbar {
                          height: 80px;
                          min-height:25px !important;
                        }
                      .navbar-nav > li > a, .navbar-brand {
                            padding-top:1px !important; 
                            padding-bottom:1px !important;
                            height: 80px;
                            }'))),
  navbarPage(
    # embed Maplet logo and title
    title = div(img(src='logo.png',
                    style="float:left; margin-top: 5px; padding-right:20px;padding-bottom:5px",
                    height = 60),
                "Krumsiek Lab",
                tags$script(HTML("var header = $('.navbar > .container-fluid');header.append('<div style=\"float:right\"><a href=\"https://weill.cornell.edu\"><img src=\"wcm2.png\" alt=\"logo\" style=\"float:right;height:50px;margin-top: 10px; padding-right:1px; \"> </a></div>');console.log(header)")),
                windowTitle = "Maplet"),
    # sticky tabs while scrolling main panel
    position = c("fixed-top"), 
    
    # Define layout of Module-Real-Time Pipeline(coded as mod6) ----------------------------------------------------
    tabPanel(HTML(paste("Real-Time", "Pipeline", sep = "<br/>")), 
             # Sidebar layout with input and output definitions ----
             dashboardPage(
               dashboardHeader(disable = TRUE),
               dashboardSidebar(disable = TRUE),
               dashboardBody(
                 sidebarLayout(
                   # Sidebar panel for inputs ----
                   sidebarPanel(
                     id = "mod6_panel1",
                     style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; width: 80%; ",
                     tags$p(
                       HTML("<b>Real-Time Pipeline Module</b> starts with original data, creates a pipeline and download it to local."
                       )),
                     tags$p(
                       HTML("Pipeline is constrained to run in order of <b>Data Loading->Preprocessing->Differential Analysis</b> and no section should be skipped."
                       )),
                     tags$p(
                       HTML("The result SE object is dependent on the <b>instant parameters</b>."
                       )),
                     # Input: Select a file ----
                     fileInput("file1", "Uploading File",
                               multiple = FALSE,
                               accept = c(".xlsx"),
                               width = "300px"),
                     # Input: Checkbox if file has header ----
                     checkboxInput("header", "Header", TRUE),
                     
                     tags$hr(),
                     
                     box(solidHeader = T, collapsible = T, collapsed = TRUE,
                         title="Sheets for Dimensions", width = "220px",
                         checkboxInput("mod6_assay_in_row", "Samples in rows?", FALSE),
                         tags$p(HTML("Assay sheet:")),
                         uiOutput("mod6_assay_sheet"),
                         tags$p(HTML("rowData sheet:")),
                         uiOutput("mod6_rowdata_sheet"),
                         tags$p(HTML("colData sheet:")),
                         uiOutput("mod6_coldata_sheet"),
                         tags$p(HTML("Click to show the original data, but this investigation is not necessary for pipeline.")),
                         actionButton("mod6_go", "Investigate")
                     ),
                     
                     tags$hr(),
                     
                     box(solidHeader = T, collapsible = T, collapsed = TRUE,
                         title="Data Loading", width = "220px",
                         tags$p(HTML("ID column in assay:")),
                         uiOutput("mod6_assay_id_column"),
                         tags$p(HTML("ID column in rowData:")),
                         uiOutput("mod6_rowdata_id_column"),
                         tags$p(HTML("ID column in colData:")),
                         uiOutput("mod6_coldata_id_column"),
                         tags$p(HTML("Run to see log text of data loading.")),
                         actionButton("mod6_go_load", "Run", width = "110px")
                     ),
                     
                     tags$hr(),
                     
                     box(solidHeader = T, collapsible = T, collapsed = TRUE,
                         title="Preprocessing", width = "220px",
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
                         tags$p(HTML("Run to see log text of data loading and preprocessing. This step may cost a few seconds to run.")),
                         actionButton("mod6_go_preprocess", "Run", width = "110px")
                     ),
                     
                     tags$hr(),
                     
                     box(solidHeader = T, collapsible = T, collapsed = TRUE,
                         title="Differential Analysis", width = "220px",
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
                         tags$p(HTML("Run to see log text of data loading, preprocessing and differential analysis. This step may cost a few seconds to run.")),
                         actionButton("mod6_go_differ", "Run", width = "110px")
                     )
                   ),
                   
                   # Main panel for displaying outputs ----
                   mainPanel(
                     id = "mod6_panel2", 
                     style = "overflow-y: auto; max-height: 100vh; position: absolute; left: 28%",
                     br(), 
                     br(), 
                     br(), 
                     # Output: Data file ----
                     tags$p(HTML("Downloading SE.Rdata may cost more than one minute. Please wait for the prompt.")),
                     downloadButton("download_se", "Download result SE .Rdata"),
                     br(),
                     br(),
                     uiOutput("mod6_main_panel")
                   )
                 )
               )
             )
    ),
    
    # Define layout of Module-Annotations Explorer(coded as mod5) ----------------------------------------------------
    tabPanel(HTML(paste("Annotations", "Explorer", sep = "<br/>")), 
             sidebarLayout(
               sidebarPanel(id = "mod5_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>Annotations Explorer Module</b> creates tables, distribution plots, or other graphics to explore the SE object."
                              )),
                            radioButtons("mod5_dimension", "Select one dimension:", 
                                         choices = list("Column Data" = "col", 
                                                        "Row Data" = "row")
                            ),
                            br(),
                            uiOutput("mod5_dimension_ui"),
                            br(),
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                              )),
                            br(),
                            # delay the output
                            actionButton("mod5_go", "Update")
               ), 
               mainPanel(id = "mod5_panel2", 
                         br(), 
                         br(), 
                         br(), 
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         uiOutput("mod5_output_ui")
               )
             )
    ),
    
    # Define layout of Module-2D Projection(coded as mod3) ----------------------------------------------------
    
    tabPanel(HTML(paste("2D", "Projection", sep = "<br/>")), 
             sidebarLayout(
               sidebarPanel(id = "mod3_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>2D Projection Module</b> generates an interactive 2D projection of PCA/UMAP."
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
    
    # Define layout of Module-All Results Explorer(coded as mod1) ----------------------------------------------------
    
    tabPanel(HTML(paste("All Results", "Explorer", sep = "<br/>")), 
             sidebarLayout(
               sidebarPanel(id = "mod1_panel1",
                            # sidebar auto-scrolling with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>All Results Explorer Module</b> extracts all the result objects one at a time."
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
                            tags$p(
                              HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection. Some plots such as box plot or multiple plots may cost dozens of seconds to show up."
                              )),
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
    
    # Define layout of Module-Feature Results Explorer(coded as mod4) --------------------------------------------------
    
    tabPanel(HTML(paste("Feature Results", "Explorer", sep = "<br/>")),
             sidebarLayout(
               sidebarPanel(
                 id = "mod4_panel1",
                 # sidebar autoscroll with main panel
                 style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                 tags$p(
                   HTML(
                     "<b>Feature Results Explorer Module</b> requires collection on all statistical results in a table given one metabolite name."
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
                 tags$p(
                   HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                   )),
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
             )
             ),
    
    # Define layout of Module-Pathway Results Explorer(coded as mod2) ----------------------------------------------------
    
    tabPanel(HTML(paste("Pathway Results", "Explorer", sep = "<br/>")),
             sidebarLayout(
               sidebarPanel(
                 id = "mod2_panel1",
                 # sidebar auto-scrolling with main panel
                 style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                 tags$p(
                   HTML(
                     "<b>Pathway Results Explorer Module</b> displays a series of interactive plots at different granularities given a SE and a statname."
                   )
                 ),
                 tags$p(
                   HTML(
                     "StatsBar plot -> Equalizer/Volcano plot -> Box/Scatter plot."
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
                 tags$hr(),
                 radioButtons(
                   "mod2.plot3",
                   "Select plot3 type:",
                   choices = list("Box" = "box",
                                  "Scatter" = "scatter"),
                   selected  = "scatter"
                 ),
                 checkboxInput("mod2.categorical", 
                               "Treat as categorical", 
                               value = FALSE
                 ),
                 br(),
                 tags$p(
                   HTML("<b>Hint:<br></b>Outputs are delayed untill you click 'UPDATE' button after selection."
                   )),
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
             ))
  )
)

################################################################################
################ Define server logic required to draw outputs ##################
################################################################################

server <- function(input, output) {
  
  # Define rendering logic of control widgets in Module-Real-Time Pipeline (coded as mod6)----------------------
  # control widget of selecting file
  output$mod6_assay_sheet <- renderUI({
    req(input$file1)
    selectInput("assay_sheet", label = NULL,
                width = "220px",
                choices = getSheetNames(as.character(input$file1$datapath))
    )
  })
  # control widget of dimensions
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
  # control widget of data loading
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
  # control widget of preprocessing
  output$mod6_pre_sample_color_column <- renderUI({
    selectInput("pre_sample_color_column", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_batch_column <- renderUI({
    selectInput("pre_batch_column", label = NULL,
                width = "220px",
                selected=NULL,
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_pca_color_column <- renderUI({
    selectInput("pre_pca_color_column", label = NULL,
                width = "220px",
                multiple=TRUE,
                selected=NULL,
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_heatmap_anno_column <- renderUI({
    selectInput("pre_heatmap_anno_column", label = NULL,
                width = "220px",
                multiple=TRUE,
                selected=NULL,
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_pre_heatmap_anno_row <- renderUI({
    selectInput("pre_heatmap_anno_row", label = NULL,
                width = "220px",
                multiple=TRUE,
                selected=NULL,
                choices = colnames(df_rowdata())
    )
  })
  # control widget of differential analysis
  output$mod6_outcome <- renderUI({
    selectInput("outcome", label = NULL,
                width = "220px",
                choices = colnames(df_coldata())
    )
  })
  
  output$mod6_group_col_barplot <- renderUI({
    selectInput("group_col_barplot", label = NULL,
                width = "220px",
                selected=NULL,
                choices = colnames(df_rowdata())
    )
  })
  
  output$mod6_color_col_barplot <- renderUI({
    selectInput("color_col_barplot", label = NULL,
                width = "220px",
                selected=NULL,
                choices = colnames(df_rowdata())
    )
  })
  
  # Define rendering logic of outputs in Module-Real-Time Pipeline(coded as mod6) ------------------------------
  # record the file path of selected file
  mod6_filepath <- 
    eventReactive(input$mod6_go, ## delayed output
                  {c(input$file1$datapath)
                  })
  
  # print table when clicking "investigate" button
  observeEvent(input$mod6_go,{
    output$mod6_main_panel <- renderUI({
      list(dataTableOutput("mod6_assay"),
           br(),br(),
           dataTableOutput("mod6_rowdata"),
           br(),br(),
           dataTableOutput("mod6_coldata"))
    })
  })
  
  # render logic of the table
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
  
  observeEvent(input$mod6_go_load,{
    # define main panel for loading section
    output$mod6_main_panel <- renderUI({
      tagAppendAttributes(verbatimTextOutput("log_load"), 
                          style="white-space:pre-wrap;")
    })
  })
  
  observeEvent(input$mod6_go_preprocess,{
    # define main panel for preprocessing section
    output$mod6_main_panel <- renderUI({
      tagAppendAttributes(verbatimTextOutput("log_preprocess"), 
                          style="white-space:pre-wrap;")
    })
  })
  
  observeEvent(input$mod6_go_differ,{
    # define main panel of differential analysis
    output$mod6_main_panel <- renderUI({
      tagAppendAttributes(verbatimTextOutput("log_differ"), 
                          style="white-space:pre-wrap;")
    })
  })
  
  # get loading SE
  D_load <- reactive({
    ## loading D
    file_data <- as.character(input$file1$datapath)
    D <-
      mt_load_xls(file=file_data, 
                  sheet=input$assay_sheet, 
                  samples_in_row=input$mod6_assay_in_row, 
                  id_col=input$assay_id_column) %>%
      mt_anno_xls(file=file_data, 
                  sheet=input$rowdata_sheet,
                  anno_type="features", 
                  anno_id_col=input$rowdata_id_column, 
                  data_id_col = "name") %>%
      mt_anno_xls(file=file_data, 
                  sheet=input$coldata_sheet, 
                  anno_type="samples", 
                  anno_id_col =input$coldata_id_column, 
                  data_id_col ="sample") %>%
      mt_reporting_data() %>%
      {.}
    ## return D
    D
  })
  
  # get proprocessing SE
  D_preprocess <- reactive({
    ## preprocessing D
    D <- D_load() %>%
      mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
      mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
      mt_plots_missingness(feat_max=input$mod6_filter_feat_max,samp_max = input$mod6_filter_sample_max) %>%
      mt_pre_filter_missingness(feat_max = input$mod6_filter_feat_max, samp_max = input$mod6_filter_sample_max) %>%
      mt_plots_missingness(feat_max=input$mod6_filter_feat_max, samp_max = input$mod6_filter_sample_max) %>%
      mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
      mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
      mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "Original", plot_logged = T) %>%
      {.}
    if(!is.null(input$pre_batch_column)){
      D %<>%
        mt_pre_batch_median(batch_col = input$pre_batch_column)
    }
    D <- D %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After batch correction", plot_logged = T) %>%
      mt_pre_norm_quot(feat_max = input$mod6_feat_max_norm) %>%
      mt_plots_dilution_factor(in_col=input$pre_sample_color_column) %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After normalization", plot_logged = T) %>%
      mt_pre_trans_log() %>%
      mt_pre_impute_knn() %>%
      mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After imputation", plot_logged = T) %>%
      mt_pre_outlier_detection_univariate() %>%
      mt_reporting_data() %>%
      mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
      {.}
    ## add PCA/UMAP plots
    lapply(input$pre_pca_color_column, function(x){
      D <<- D %>%
        mt_plots_pca(scale_data = T, title = sprintf("scaled PCA - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
        mt_plots_umap(scale_data = T, title = sprintf("scaled UMAP - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
        {.}
    }) %>% invisible
    ## add heatmap
    D %<>%
      mt_plots_heatmap(scale_data = T, annotation_col = input$pre_heatmap_anno_column, annotation_row = input$pre_heatmap_anno_row,
                       clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
      {.}
    ## return D
    D
  })
  
  ## get differential analysis SE
  D_differ <- reactive({
    # Differential analysis D
    D <- D_preprocess() %>%
      mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
      diff_analysis_func(var=input$outcome,
                         binary=input$mod6_outcome_binary,
                         analysis_type=input$mod6_analysis_type, 
                         mult_test_method=input$mod6_mult_test_method,
                         alpha=input$mod6_sig_threshold,
                         group_col_barplot=input$group_col_barplot,
                         color_col_barplot=input$color_col_barplot) %>%
      {.}
    ## return D
    D
  })
  
  # render logic of the log text of data loading
  output$log_load <- renderPrint({
    get_log_text(D_load())
  })
  # render logic of the log text of preprocessing
  output$log_preprocess <- renderPrint({
    # loading log
    text_load <- get_log_text(D_load())
    # preprocessing log
    text_preprocess <- get_log_text(D_preprocess())
    # paste log text
    str <- paste(text_load, text_preprocess, sep = "\n")
    cat(str)
  })
  
  # render logic of the log text of differential analysis
  output$log_differ <- renderPrint({
    # loading log
    text_load <- get_log_text(D_load())
    # preprocessing log
    text_preprocess <- get_log_text(D_preprocess())
    # differential analysis log
    text_differ <- get_log_text(D_differ())
    # paste log text
    str <- paste(text_load, text_preprocess, text_differ, sep = "\n")
    cat(str)
  })
  
  # download SE button
  # https://mastering-shiny.org/action-transfer.html
  output$download_se <- downloadHandler(
    filename = function() {
      paste0("SE_", Sys.Date(), ".Rdata")
    },
    content = function(fname) {
      ## loading D
      file_data <- as.character(input$file1$datapath)
      D <-
        mt_load_xls(file=file_data, 
                    sheet=input$assay_sheet, 
                    samples_in_row=input$mod6_assay_in_row, 
                    id_col=input$assay_id_column) %>%
        mt_anno_xls(file=file_data, 
                    sheet=input$rowdata_sheet,
                    anno_type="features", 
                    anno_id_col=input$rowdata_id_column, 
                    data_id_col = "name") %>%
        mt_anno_xls(file=file_data, 
                    sheet=input$coldata_sheet, 
                    anno_type="samples", 
                    anno_id_col =input$coldata_id_column, 
                    data_id_col ="sample") %>%
        mt_reporting_data() %>%
        {.}
      ## preprocessing D
      D <- D %>%
        mt_reporting_heading(heading = "Preprocessing", lvl=1) %>%
        mt_reporting_heading(heading = "Filtering", lvl = 2) %>%
        mt_plots_missingness(feat_max=input$mod6_filter_feat_max,samp_max = input$mod6_filter_sample_max) %>%
        mt_pre_filter_missingness(feat_max = input$mod6_filter_feat_max, samp_max = input$mod6_filter_sample_max) %>%
        mt_plots_missingness(feat_max=input$mod6_filter_feat_max, samp_max = input$mod6_filter_sample_max) %>%
        mt_anno_missingness(anno_type = "samples", out_col = "missing") %>%
        mt_anno_missingness(anno_type = "features", out_col = "missing") %>%
        mt_reporting_heading(heading = "Normalization", lvl = 2) %>%
        mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "Original", plot_logged = T) %>%
        {.}
      if(!is.null(input$pre_batch_column)){
        D %<>%
          mt_pre_batch_median(batch_col = input$pre_batch_column)
      }
      D <- D %>%
        mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After batch correction", plot_logged = T) %>%
        mt_pre_norm_quot(feat_max = input$mod6_feat_max_norm) %>%
        mt_plots_dilution_factor(in_col=input$pre_sample_color_column) %>%
        mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After normalization", plot_logged = T) %>%
        mt_pre_trans_log() %>%
        mt_pre_impute_knn() %>%
        mt_plots_sample_boxplot(color=!!sym(input$pre_sample_color_column), title = "After imputation", plot_logged = T) %>%
        mt_pre_outlier_detection_univariate() %>%
        mt_reporting_data() %>%
        mt_reporting_heading(heading = "Global Statistics", lvl = 1) %>%
        {.}
      ## add PCA/UMAP plots
      lapply(input$pre_pca_color_column, function(x){
        D <<- D %>%
          mt_plots_pca(scale_data = T, title = sprintf("scaled PCA - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
          mt_plots_umap(scale_data = T, title = sprintf("scaled UMAP - %s",x), color=!!sym(x), size=2.5, ggadd=scale_size_identity()) %>%
          {.}
      }) %>% invisible
      ## add heatmap
      D %<>%
        mt_plots_heatmap(scale_data = T, annotation_col = input$pre_heatmap_anno_column, annotation_row = input$pre_heatmap_anno_row,
                         clustering_method = "ward.D2", fontsize = 5, cutree_rows = 3, cutree_cols = 3, color=gplots::bluered(101)) %>%
        {.}
      # Differential analysis D
      D <- D %>%
        mt_reporting_heading(heading = "Statistical Analysis", lvl = 1) %>%
        diff_analysis_func(var=input$outcome,
                           binary=input$mod6_outcome_binary,
                           analysis_type=input$mod6_analysis_type, 
                           mult_test_method=input$mod6_mult_test_method,
                           alpha=input$mod6_sig_threshold,
                           group_col_barplot=input$group_col_barplot,
                           color_col_barplot=input$color_col_barplot) %>%
        {.}
      # write Rdata to local
      save(D, file=fname)
    }
  )
  
  # Define rendering logic of control widgets in Module-Annotations Explorer(coded as mod5) ----------------------
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
                      tags$hr(),
                      selectInput("mod5_var2_select", 
                                  "Select the secondary variable:", 
                                  choices = names(colData(D)),
                                  selected = "sample",
                                  width = "220px"),
                      checkboxInput("mod5_var2_type", 
                                    "Continuous", 
                                    value = TRUE),
                      tags$hr(),
                      selectInput("mod5_select_hover", 
                                  "Select hovering text:", 
                                  choices = names(colData(D)),
                                  selected = names(colData(D))[1],
                                  width = "220px",
                                  multiple=TRUE)
           ),
           "row"=selectInput("mod5_rowdata_plot", 
                             "Select one plot for row data:", 
                             choices = c("SUPER_PATHWAY"),
                             width = "220px")
    )
  })
  
  output$mod5_output_ui <- renderUI({
    switch(input$mod5_dimension,
           "col"=list(downloadButton("mod5_download_plotly", "download plotly"),
                      plotlyOutput('mod5_plot', height = 600)),
           "row"=list(fluidRow(
             splitLayout(style = "border: 1px", cellWidths = c(1000, 1000), 
                         downloadButton("mod5_download_plotly", "download plotly"), 
                         downloadButton("mod5_download_plotly2", "download plotly")
             )
           ),
           fluidRow(
             splitLayout(style = "height:600px; border: 1px", cellWidths = c(1000, 1000), 
                         plotlyOutput('mod5_plot', height = 600), 
                         plotlyOutput('mod5_plot2', height = 600)
             )
           ))
    )
  })
  
  # Define rendering logic of outputs in Module-Annotations Explorer(coded as mod5) ------------------------------
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
                                            mod5_scatter(D, x=mod5_input()[3], 
                                                         y=mod5_input()[1], 
                                                         hover = input$mod5_select_hover)
                                          } else if(mod5_input()[2]==TRUE & mod5_input()[4]==FALSE) {
                                            mod5_boxplot(D, x=mod5_input()[3], 
                                                         x_cate = FALSE,
                                                         y=mod5_input()[1],
                                                         y_cate = TRUE,
                                                         fill=mod5_input()[3], 
                                                         hover=input$mod5_select_hover)
                                          } else if(mod5_input()[2]==FALSE & mod5_input()[4]==TRUE) {
                                            mod5_boxplot(D, x=mod5_input()[1], 
                                                         x_cate = FALSE,
                                                         y=mod5_input()[3],
                                                         y_cate = TRUE,
                                                         fill=mod5_input()[1], 
                                                         hover=input$mod5_select_hover)
                                          } else {
                                            mod5_barplot(D, x=mod5_input()[3], 
                                                         fill=mod5_input()[1], 
                                                         hover = input$mod5_select_hover)
                                          },
                                        "row"=
                                          rowData(D) %>%
                                          data.frame %>%
                                          dplyr::rename(var=mod5_input()[5]) %>%
                                          dplyr::group_by(var) %>%
                                          dplyr::summarise(count=n()) %>%
                                          plot_ly(labels = ~var, 
                                                  values = ~count, 
                                                  type = 'pie',
                                                  textposition = 'inside',
                                                  source="mod5-click",
                                                  title="<b>Distribution of Super Pathway</b>") %>% 
                                          layout(autosize = F, width = 1000, height = 500,
                                                 uniformtext=list(minsize=12, mode='hide'),
                                                 legend = list(x = 1,
                                                               y = .5,
                                                               tracegroupgap = 5)
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
  
  ## to see the stored data of clicking
  # output$info <- renderPrint({
  #   d5 <- event_data("plotly_click", source = "mod5-click")
  #   if(!is.null(d5)){
  #     d5
  #   }
  # })
  
  output$mod5_plot2 <- renderPlotly({
    d5 <- event_data("plotly_click", source = "mod5-click")
    pie_dat <- as.data.frame(rowData(D))
    
    if (!is.null(d5)){
      lvls <- rev(pie_dat$SUPER_PATHWAY)
      label <- lvls[round(as.numeric(d5$pointNumber))+1]
      
      session_store$mod5_plot2 <- 
        pie_dat[pie_dat$SUPER_PATHWAY == label, ] %>%
        dplyr::rename(var="SUB_PATHWAY") %>%
        dplyr::group_by(var) %>%
        dplyr::summarise(count=n()) %>%
        plot_ly(labels = ~var, 
                values = ~count, 
                type = 'pie', 
                textposition = 'inside',
                title=paste0("<b>Distribution of Sub Pathway in Specified Super Pathway - </b>", label)
        ) %>%
        layout(autosize = F, width = 1000, height = 500,
               uniformtext=list(minsize=12, mode='hide'),
               legend = list(x = 1,
                             y = .5,
                             tracegroupgap = 5)
        )
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
  
  # Define rendering logic of control widgets in Module-Annotations Explorer(coded as mod3) ------------------------
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
  
  # Define rendering logic of outputs in Module-Annotations Explorer(coded as mod3) --------------------------------
  
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
  
  # Define rendering logic of control widgets in Module-All Results Explorer(coded as mod1) ------------------------
  
  # create stat_name list dependent on radio button
  output$mod1_select_statname_ui <- renderUI({
    selectInput("mod1_select_statname", "Select one stat name:",
                width = "220px",
                choices = dplyr::distinct(obj_name[obj_name$V1==input$mod1_radio, ], stat_name)$stat_name
    )
  })
  
  # create object list dependent on radio button and stat_name
  output$mod1_select_object_ui <- renderUI({
    if (input$mod1_radio=="stats"){
      NULL
    } else {
      selectInput("mod1_select_object", "Select one object:",
                  width = "220px",
                  choices = dplyr::distinct(obj_name[obj_name$stat_name==input$mod1_select_statname&obj_name$V1==input$mod1_radio, ], V2)$V2
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
    # assign a data frame of all the object names of box plots
    # filter() cannot run in Shiny, use subset() instead
    box_obj_name <- subset(obj_name, V1=="plots"&V2=="box")
    box_output_order <- box_obj_name %>%
      dplyr::mutate(order=seq(from=1, to=n()))
    
    if(input$mod1_select_statname %in% box_output_order$stat_name){
      box_output_order[box_output_order$stat_name==input$mod1_select_statname, ]$order
    } else {
      1
    }
  })
  
  # create reactive inputs list
  mod1_input_object <- eventReactive(input$mod1_go, ## delayed output
                                     {c(input$mod1_radio,
                                        input$mod1_select_statname,
                                        input$mod1_select_object)}
  )
  
  # Define rendering logic of outputs in Module-All Results Explorer(coded as mod1) --------------------------------
  
  # Insert the right number of plot output objects into UI
  output$mod1_output_plot <- renderUI({
    ## limit plots to specified stat_name
    obj_name <- subset(obj_name, V1==mod1_input_object()[1])
    obj_name <- subset(obj_name, V2==mod1_input_object()[3])
    output_order <- obj_name %>%
      dplyr::mutate(order=seq(from=1, to=n()))
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
        as.numeric(plots[[row_n]]$output2)*150
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
    # get the max number of objects
  num_df <- subset(obj_name, V1=="plots")
  num_df <- num_df %>%
    dplyr::group_by(V2, stat_name) %>%
    dplyr::summarise(cnt_sum=sum(cnt))
  max_plot <- max(num_df$cnt_sum)
  
  for (i in 1:max_plot) {
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
          dplyr::mutate(order=seq(from=1, to=n()))
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
      dplyr::select(c(2, 20:26))
    
    ## put interested columns ahead
    table <- if ('term' %in% names(table)) {
      table %>%
        dplyr::select(name, statistic, p.value, p.adj, term, dplyr::everything()) %>%
        ## scientific notation
        dplyr::mutate(statistic=formatC(statistic, format = "E", digits = 2),
                      p.value=formatC(p.value, format = "E", digits = 2),
                      p.adj=formatC(p.adj, format = "E", digits = 2),
                      estimate=formatC(estimate, format = "E", digits = 2),
                      std.error=formatC(std.error, format = "E", digits = 2)
        )
    } else {
      table %>%
        dplyr::select(name, statistic, p.value, p.adj, dplyr::everything())
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
  
  # Define rendering logic of outputs in Module-Feature Results Explorer(coded as mod4) --------------------------------
  
  # Module 4: general reactive stats table
  mod4_metabolite_table <-
    eventReactive(input$mod4_go,
                  {
                    table <- data.frame()
                    # assign an object of all stats tables
                    table_stats <- mtm_res_get_entries(D, c("stats", "univ"))
                    # assign an object of all stats plots
                    plot_stats <- mtm_res_get_entries(D, c("plots", "stats"))
                    for (i in 2:length(table_stats)) {
                      tab <- table_stats[[i]]$output$table %>%
                        dplyr::mutate(`stat name` = plot_stats[[i - 1]]$args$stat_list)
                      table <- rbind(table, tab)
                    }
                    table <- table %>%
                      dplyr::select(var, statistic, p.value, p.adj, `stat name`, estimate, std.error) %>%
                      dplyr::mutate(
                        statistic = formatC(statistic, format = "E", digits = 2),
                        p.value = formatC(p.value, format = "E", digits = 2),
                        p.adj = formatC(p.adj, format = "E", digits = 2),
                        estimate = formatC(estimate, format = "E", digits = 2),
                        std.error = formatC(std.error, format = "E", digits = 2)
                      ) %>%
                      dplyr::filter(var == input$mod4_metabolite) %>%
                      dplyr::rename("name" = var)
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
  
  observe({
    if (!is.null(input$mod4_table_rows_selected)) {
      session_store$mod4.tb.row <- input$mod4_table_rows_selected
    }
  })
  
  # mod4: extract the stat_name
  stat_name_selected <- reactive({
    mod4_metabolite_table() %>% 
      dplyr::slice(round(as.numeric(session_store$mod4.tb.row))) %>%
      dplyr::select(`stat name`)
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
    highlight_point <- data_vol[data_vol$isSelected==TRUE, ]
    
    plot <- data_vol %>%
      ggplot(aes(x = statistic, y = p.value, color = isSelected, label = name)) +
      geom_point() +
      geom_point(data=highlight_point, size = 3) +
      scale_y_continuous(trans = reverselog_trans(10),
                         breaks = scales::trans_breaks("log10", function(x) 10^x),
                         labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      labs(y = "p-value (10^(-y))") +
      ggtitle(paste0(stat_name_selected(), "-", isSelected)) +
      scale_color_manual(values=c("#999999", "red"))
    
    session_store$mod4.vol <- ggplotly(plot, source = "mod4_sub_vol") %>%
      layout(legend = list(orientation = 'h',
                           xanchor = "center",
                           x = 0.5,
                           y = -0.2,
                           title = list(text='<b> isSelected </b>')))
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
      
      # set the column curveNumber by color legend
      isSelected <- input$mod4_metabolite
      data_vol[, "curveNumber"] <- ifelse(data_vol$var==isSelected, 1, 0)
      data_vol_true <- data_vol[data_vol$curveNumber==1, ]
      data_vol_false <- data_vol[data_vol$curveNumber==0, ]
      
      # By using click info (curveNumber & ponitNumber) to get the metabolite name
      metabolite <- ifelse(d$curveNumber == 1,
                           data_vol_true[d$pointNumber + 1, ]$var[1],
                           data_vol_false[d$pointNumber + 1, ]$var[1])
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
          ggplot(aes(x = !!sym(term), y = value)) +
          geom_point(size = 3) +
          geom_smooth(method = "lm", se = T, color = "black") + 
          ggtitle(metabolite)
      } else {
        plot <- data %>%
          ggplot(aes(x = !!sym(term), y = value)) +
          geom_boxplot() +
          geom_jitter(size = 3, width = 0.2) +
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
  
  # Define rendering logic of outputs in Module-Pathway Results Explorer(coded as mod2) --------------------------------
  # Module 2: create reactive inputs list
  mod2_input_object <- eventReactive(input$mod2_go, 
                                     {c(input$mod2.stat,
                                        input$mod2.plot1,
                                        input$mod2.plot2,
                                        input$mod2.plot3)}
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
    
    # equalizer/bar -> bar/null -> plot
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
    d <- event_data("plotly_click", source = "sub_bar")
    legend_name <- paste0("p.adj < ", alpha)
    
    if (!is.null(d)) {
      # D:SE object, inputs: sidebar value, legend_name: legend name
      # d: click info for bar plot, pwvar: SUB_PATWAY/PATTHWAY, alpha: significant value (ex. p.adj < 0.1)
      plot <- mod2_plot_vol(D, inputs, legend_name, d, pwvar, alpha)
    } else {
      plot <- mod2_plot_vol(D, inputs, legend_name, NULL, pwvar, alpha)
    }
    
    session_store$mod2.vol <- ggplotly(plot, source = "sub_vol") %>%
      layout(legend = list(orientation = 'h',
                           xanchor = "center",
                           x = 0.5,
                           y = -0.2,
                           title = list(text=paste0('<b> ', legend_name, ' </b>'))))
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
    d <- event_data("plotly_click", source = "sub_bar")
    
    if (inputs[2] == "null") {
      # D:SE object, inputs: sidebar value, rd: pathway annotations
      # alpha: significant value (ex. p.adj < 0.1), pwvar: SUB_PATWAY/PATTHWAY,
      # path_name: pathway name for equalizer plot, d: click info for bar plot
      plot <- mod2_plot_eq(D, inputs, rd, alpha, pwvar, input$mod2.equal.path, NULL)
    } else {
      plot <- mod2_plot_eq(D, inputs, rd, alpha, pwvar, NULL, d)
    }
    
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
    
    d.bar <- event_data("plotly_click", source = "sub_bar")
    d.eq <- event_data("plotly_click", source = "sub_eq")
    d.vol <- event_data("plotly_click", source = "sub_vol")
    
    plot <- mod2_plot_box_scatter(D, # SE object
                                  inputs, # sidebar inputs
                                  d.bar, # click info for bar plot
                                  d.eq, # click info for equalizer plot
                                  d.vol, # click info for volcano plot
                                  rd, # pathway annotations
                                  pwvar, # pathway annotation column 
                                  input$mod2.equal.path, # pathway name if plot2 is "equalizer"
                                  alpha, # significant value (ex. p.adj < 0.1)
                                  input$mod2.categorical, # if treated categorical
                                  data) # data for box/scatter plot
    
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
  
}

# Run the application 
shinyApp(ui = ui, server = server)