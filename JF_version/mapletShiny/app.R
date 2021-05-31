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

# load packages
library(shiny)
library(shinyWidgets)
library(maplet)
# devtools::install_github(repo="krumsieklab/maplet@v1.0.1", subdir="maplet")
library(tidyverse)
library(DT)
library(plotly)

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
    # customize hover text
    # hover_text <- paste0(names(data.frame(colData(D)))[as.numeric(hover)], ": ", 
    #                     data.frame(colData(D))[[as.numeric(hover)]])
    # draw ggplot
    p <- ggplot(data = df, 
                aes_string(
                  x = sprintf("PC%d", pc1), 
                  y = sprintf("PC%d", pc2),
                  color = color,
                  text=hover)
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
    colnames(df)[1:2] <- c(sprintf("PC%d", pc1), 
                           sprintf("PC%d", pc2)
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
                aes_string(
                  x = sprintf("PC%d", pc1), 
                  y = sprintf("PC%d", pc2),
                  color=color,
                  text=hover)
                ) + 
      geom_point() + 
      xlab(pc1name) + 
      ylab(pc2name) + 
      ggtitle(plot_title)
  }
  
  # draw plotly
  ggplotly(p, tooltip = "text") %>%
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
  # customize hover text
  # hover_text <- paste0(names(data.frame(colData(D)))[as.numeric(hover)], ": ", 
  #                      data.frame(colData(D))[[as.numeric(hover)]])
  # draw ggplot
  p <- ggplot(data = df,
              aes_string(x = "comp1", 
                         y = "comp2",
                         color=color,
                         text=hover)
              ) + 
    geom_point() + 
    xlab("comp 1") + 
    ylab("comp 2") + 
    ggtitle(plot_title) +
    ## reactive legend title
    labs(color = color)
  
  # draw plotly
  ggplotly(p, tooltip = "text") %>% 
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = .5,
                         y = -.2,
                         tracegroupgap = 5),
           autosize = TRUE
    )
  
}

# define equilizer plot object

# define pathway annotation column (extracted from corresponding stat_bar plot)
pwvar <- "pathway"
# define threshold for significance (extracted from corresponding stat_bar plot)
alpha <- 0.2
# define stat_names (these should be extracted dynamically from the SE object
comp <- stat_name %>% filter(!is.na(stat_name)) %>% unique() %>% unlist()

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
# generate equalizer plots for all stat_names and all pathways
eqplot <- lapply(comp %>% {names(.)=.;.}, function(comparison){
  # add pathway annotations to results
  res <- maplet::mtm_get_stat_by_name(D, comparison) %>%
    dplyr::left_join(rd, by=c("var"="name")) %>%
    dplyr::mutate(x = sign(statistic)*log10(p.adj)) %>%
    dplyr::filter(!is.na(BIOCHEMICAL))
  # colors
  clrs <- c("#9494FF","red")
  # compute multiple testing correction line
  sel <- res %>% 
    dplyr::filter(p.adj < alpha) 
  if(nrow(sel)>0){
    xfine <- sel %>% .$p.adj %>% max(., na.rm = T)
  } else {
    xfine <- Inf
  }
  # get pathway names
  pw <- res %>%
    dplyr::pull(!!dplyr::sym(pwvar)) %>% 
    unique
  # create equalizer plots
  lapply(pw %>% {names(.)=.;.}, function(y){
    # select metabolites in pw x
    df <- res %>%
      dplyr::filter(!!sym(pwvar)==y)
    # x axis limits
    a = max(abs(df$x))+0.3
    ggplot(df, aes(x = x, y = BIOCHEMICAL)) +
      geom_vline(xintercept = 0, color ="gray") +
      (if(!is.infinite(xfine)){geom_vline(xintercept = c(-log10(xfine),log10(xfine)), color="red", alpha=0.4)}) +
      (if(!is.infinite(xfine)){ggtitle(sprintf("Differential Metabolites at alpha %.2f",alpha))}else{ggtitle(sprintf("No significant results at alpha %.2f", alpha))}) +
      geom_point(pch = 22, fill = clrs[1], size = 3) +
      facet_grid(as.formula(sprintf("%s~.",pwvar)), scales = "free_y", space = "free_y") +
      theme(strip.background =element_rect(fill=NA),
            strip.text = element_text(colour = 'black', face = "bold"),
            strip.text.y = element_text(angle = 0, hjust = 0),
            panel.grid.major.y = element_line(color ="gray"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank(),
            panel.background = element_rect(fill=NA, color ="black")) +
      ylab("") +
      xlab("sign(statistic)*log10(p.adj)") +
      scale_x_continuous(limits = c(-a,a))
  })
})


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
                   "mod2.plot2",
                   "Select plot2 type:",
                   choices = list("Equalizer" = "equalizer",
                                  "Volcano" = "volcano"),
                   selected  = "volcano"
                 ),
                 br(),
                 radioButtons("mod2.plot3",
                              "Select plot3 type:",
                              choices = list("Box" = "box",
                                             "Scatter" = "scatter"),
                              selected  = "scatter"
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
                 style = "overflow-y: auto; position: absolute; left: 25%",
                 downloadButton("mod2_download_plotly_bar", "download bar plot"),
                 plotlyOutput("mod2.bar", height = 600),
                 br(),
                 downloadButton("mod2_download_plotly_volcano", "download volcano plot"),
                 plotlyOutput("mod2.vol", height = 600),
                 br(),
                 downloadButton("mod2_download_plotly_scatter", "download scatter plot"),
                 plotlyOutput("mod2.box", height = 400),
                 br()
               )
             ) 
    ),

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
                         plotlyOutput('mod3_plot', height = 730)
               )
             )
    ),

  # Define layout of Module 4 ----------------------------------------------------

    tabPanel("Module 4", 
             sidebarLayout(
               sidebarPanel(id = "mod4_panel1",
                            # sidebar autoscroll with main panel
                            style = "margin-left: -25px; margin-top: 45px; margin-bottom: 5px; position:fixed; width: 20%; height: 100%;",
                            tags$p(
                              HTML("<b>Module 4</b> requires collection on all statistical results in a table given one metabolite name."
                              )),
                            tags$p(
                              HTML("When clicking on one row, it should display interactive plots following the same orders in Module 2."
                              )),
                            # select one metabolite
                            selectInput("mod4_metabolite", "Select one metabolite:",
                                        width = "220px",
                                        choices = arrange(mtm_res_get_entries(D, c("stats", "univ"))[[1]]$output$table, var)$var,
                                        selected = "pantoate"
                            ),
                            br(),
                            # delay the output
                            actionButton("mod4_go", "Update")
               ), 
               mainPanel(id = "mod4_panel2", 
                         br(), 
                         br(), 
                         br(), 
                         style = "overflow-y: auto; position: absolute; left: 25%",
                         # stats table
                         dataTableOutput('mod4_table'),
                         br(),
                         # statsBar plotly
                         downloadButton("mod4_download_plotly_bar", "download bar plot"),
                         plotlyOutput('mod4_stats_bar', height = 730),
                         br(), 
                         # volcano plotly
                         downloadButton("mod4_download_plotly_volcano", "download volcano plot"),
                         plotlyOutput('mod4_volcano', height = 730),
                         br(), 
                         # box/scatter plotly
                         downloadButton("mod4_download_plotly_scatter", "download scatter plot"),
                         plotlyOutput('mod4_box_scatter', height = 730)
               )
             )
    ),
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
                         plotlyOutput('mod5_plot', height = 600)
               )
             )
    ),
    tabPanel("Module 6", "contents")
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
    table <- data.frame(metabolite=row.names(rowData(D)), rowData(D)) %>%
      left_join(mtm_get_stat_by_name(D, mod1_input_object()[2]), 
                by=c("metabolite"="var")
                ) %>%
      select(c(4, 20:26)) %>%
      ## scientific notation
      mutate(statistic=formatC(statistic, format = "E", digits = 2),
             p.value=formatC(p.value, format = "E", digits = 2),
             p.adj=formatC(p.adj, format = "E", digits = 2)
      )
    
    ## put interested columns ahead
    table <- if ('term' %in% names(table)) {
      table %>%
        select(BIOCHEMICAL, statistic, p.value, p.adj, term, everything())
    } else {
      table %>%
        select(BIOCHEMICAL, statistic, p.value, p.adj, everything())
    }
    datatable(table,
              options = list(
                # limit number of rows
                pageLength =  10,
                lengthMenu = c(10, 20, 50),
                ## set column width
                autoWidth = TRUE,
                columnDefs = list(list(width = '50px', targets = c(2:4))),
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
  session_store <- reactiveValues()
  
  # create reactive inputs list
  mod2_input_object <- eventReactive(input$mod2_go, 
                                     {c(input$mod2.stat,
                                        input$mod2.plot2,
                                        input$mod2.plot3)}
  )
  
  # barplot
  output$mod2.bar <- renderPlotly({
    inputs <- mod2_input_object()
    plots <- mtm_res_get_entries(D, c("plots", "stats"))
    for (i in seq_along(plots)) {
      if (plots[[i]]$args$stat_list == inputs[1]) {
        plot <- plots[[i]]$output[[1]]
      }
    }
    session_store$mod2.bar <- ggplotly(plot, source = "mod2_sub_bar") %>%
      layout(dragmode = "lasso",
             legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3))
    session_store$mod2.bar
  })
  
  # download button 
  output$mod2_download_plotly_bar <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # export plotly html widget as a temp file to download.
      saveWidget(as_widget(session_store$mod2.bar), file, selfcontained = TRUE)
    }
  )
  
  # volcano plot
  
  output$mod2.vol <- renderPlotly({
    inputs <- mod2_input_object()
    d <- event_data("plotly_click", source = "mod2_sub_bar")
    if (!is.null(d)) {
      plots_bar <- mtm_res_get_entries(D, c("plots", "stats"))
      for (i in seq_along(plots_bar)) {
        if (plots_bar[[i]]$args$stat_list == inputs[1]) {
          data_bar <- plots_bar[[i]]$output[[1]]$data
        }
      }
      lvls <- rev(levels(data_bar$label))
      label <- lvls[round(as.numeric(d$y))]
      name <- data_bar[data_bar$label == label, ]$name
      
      plots_vol <- mtm_res_get_entries(D, c("plots", "volcano"))
      for (i in seq_along(plots_vol)) {
        if (plots_vol[[i]]$args$stat_name == inputs[1]) {
          data_vol <- plots_vol[[i]]$output[[1]]$data
        }
      }
      
      row_data <- rowData(D) %>% data.frame()
      names <- unlist(row_data[row_data$SUB_PATHWAY == name,]$name)
      data_vol$isSelected <- data_vol$name %in% names
      data_vol$text <- ifelse(data_vol$name %in% names, data_vol$name, "")
      
      t <- list(family = "sans serif", size = 14, color = toRGB("grey50"))
      
      plot <- data_vol %>%
        ggplot(aes(x = statistic, y = p.value, color = isSelected)) +
        geom_point(aes(text = name)) +
        scale_y_continuous(trans = reverselog_trans(10),
                           breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(y = "p-value") +
        ggtitle(paste0(inputs[1], "-", name)) +
        ggrepel::geom_text_repel(aes(label = name), max.overlaps = Inf)
      
      session_store$mod2.vol <- ggplotly(plot, source = "mod2_sub_vol") %>%
        layout(dragmode = "lasso",
               legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3)) %>%
        add_text(text=~data_vol$text,
                 textposition="top right",
                 showlegend = T)
      session_store$mod2.vol
    }
  })
  
  # download button
  output$mod2_download_plotly_volcano <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod2.vol), file, selfcontained = TRUE)
    }
  )
  
  # box-scatter plot
  output$mod2.box <- renderPlotly({
    inputs <- mod2_input_object()
    d1 <- event_data("plotly_click", source = "mod2_sub_bar")
    d2 <- event_data("plotly_click", source = "mod2_sub_vol")
    if (!is.null(d1) & !is.null(d2)) {
      plots_bar <- mtm_res_get_entries(D, c("plots", "stats"))
      for (i in seq_along(plots_bar)) {
        if (plots_bar[[i]]$args$stat_list == inputs[1]) {
          data_bar <- plots_bar[[i]]$output[[1]]$data
        }
      }
      lvls <- rev(levels(data_bar$label))
      label <- lvls[round(as.numeric(d1$y))]
      name <- data_bar[data_bar$label == label, ]$name
      
      plots_vol <- mtm_res_get_entries(D, c("plots", "volcano"))
      for (i in seq_along(plots_vol)) {
        if (plots_vol[[i]]$args$stat_name == inputs[1]) {
          data_vol <- plots_vol[[i]]$output[[1]]$data
        }
      }
      
      row_data <- rowData(D) %>% data.frame()
      names <- unlist(row_data[row_data$SUB_PATHWAY == name,]$name)
      data_vol <- data_vol[data_vol$name %in% names, ]
      name2 <- data_vol[as.numeric(d2$pointNumber) + 1, ]$name[1]
      
      rd <- rowData(D) %>%
        as.data.frame() %>%
        dplyr::mutate(var = rownames(D))
      
      stat <- maplet::mtm_get_stat_by_name(D, inputs[1]) %>%
        dplyr::inner_join(rd, by = "var")
      
      data_box <- D %>%
        maplet:::mti_format_se_samplewise() %>%
        tidyr::gather(var, value, dplyr::one_of(rownames(D)))%>%
        dplyr::rename("Age Regression" = "Age",
                      "Comparison Outcome1" = "outcome1",
                      "Comparison Outcome2" = "outcome2",
                      "Diagnosis Analysis" = "Diagnosis") %>%
        dplyr::select(dplyr::one_of(c("var","value", inputs[1]))) %>%
        dplyr::inner_join(stat[,dplyr::intersect(colnames(stat),
                                                 c('var','statistic','p.value','p.adj','name'))], by = "var") %>%
        dplyr::select(-var)
      
      data_box <- data_box[data_box$name==name2, ]
      p.value <- signif(mean(data_box$p.value), 3)
      p.adj <- signif(mean(data_box$p.adj), 3)
      plot <- data_box %>%
        ggplot(aes(x = !!sym(input$mod2.stat),
                   y = value)) +
        geom_point(aes(text = paste0("Name:", name))) +
        geom_smooth(method = "lm", se = FALSE) +
        ggtitle(paste0(input$mod2.stat, "-", name, "-", name2))
      
      session_store$mod2.box <- ggplotly(plot) %>%
        layout(dragmode = "lasso") %>%
        add_annotations(
          x = min(data_box$Age) + 5,
          y = max(data_box$value),
          text = paste0("P-value: ", p.value),
          showarrow = F
        ) %>%
        add_annotations(
          x = min(data_box$Age) + 5,
          y = max(data_box$value) - 1,
          text = paste0("P.adj: ", p.adj),
          showarrow = F
        )
      session_store$mod2.box
    }
  })
  
  # download button
  output$mod2_download_plotly_scatter <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod2.box), file, selfcontained = TRUE)
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
                  multiple=FALSE
                  )
      ),
      "pca-loadings"=list(
      checkboxInput("mod3_scale_data", "Scaled data", 
                    value = TRUE
                    ),
      selectInput("mod3_select_colData", 
                  "Select one coloring variable:", 
                  choices = names(rowData(D)),
                  selected = "BIOCHEMICAL",
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
                  selected = "BIOCHEMICAL",
                  width = "220px",
                  multiple=FALSE
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
                  multiple=FALSE
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
  
  # create reactive stats table
  mod4_metabolite_table <- 
    eventReactive(input$mod4_go,
                  {
                    table <- data.frame()
                    for (i in 2:length(table_stats)){
                      tab <- table_stats[[i]]$output$table %>%
                        mutate(`stat name`=plot_stats[[i-1]]$args$stat_list)
                      table <- rbind(table, tab)
                    }
                    table <- table %>%
                      select(var, statistic, p.value, p.adj, `stat name`, estimate, std.error) %>%
                      mutate(statistic=formatC(statistic, format = "E", digits = 2),
                             p.value=formatC(p.value, format = "E", digits = 2),
                             p.adj=formatC(p.adj, format = "E", digits = 2),
                             estimate=formatC(estimate, format = "E", digits = 2),
                             std.error=formatC(std.error, format = "E", digits = 2)
                      ) %>%
                      filter(var==input$mod4_metabolite)
                  }
    )
  
  # render stats table of Mod4
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
  
  # catch the selected row
  session_store$s <- NA
  observe({
    if(!is.null(input$mod4_table_rows_selected)){
      session_store$s <- input$mod4_table_rows_selected
    }
  })
  
  ## extract the stat_name
  stat_name_selected <- reactive({
    mod4_metabolite_table() %>% 
      slice(round(as.numeric(session_store$s))) %>%
      select(`stat name`)
  })
  
  # render statsBar plot in Mod4
  output$mod4_stats_bar <- renderPlotly({
    for (i in 1:(length(plot_stats)-1)) {
      if (plot_stats[[i]]$args$stat_list == stat_name_selected()) {
        plot <- plot_stats[[i]]$output[[1]]
      }
    }
    session_store$mod4_stats_bar <- ggplotly(plot, source = "mod4_sub_bar") %>% 
      layout(dragmode = "lasso")
    session_store$mod4_stats_bar
  })
  # download button
  output$mod4_download_plotly_bar <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      # export plotly html widget as a temp file to download.
      saveWidget(as_widget(session_store$mod4_stats_bar), file, selfcontained = TRUE)
    }
  )
  
  # render volcano plot in Mod4
  output$mod4_volcano <- renderPlotly({
    d <- event_data("plotly_click", source = "mod4_sub_bar")
    if (!is.null(d)) {
      plots_bar <- mtm_res_get_entries(D, c("plots", "stats"))
      for (i in 1:(length(plots_bar)-1)) {
        if (plots_bar[[i]]$args$stat_list == stat_name_selected()) {
          data_bar <- plots_bar[[i]]$output[[1]]$data
        }
      }
      lvls <- rev(levels(data_bar$label))
      label <- lvls[round(as.numeric(d$y))]
      name <- data_bar[data_bar$label == label, ]$name
      
      plots_vol <- mtm_res_get_entries(D, c("plots", "volcano"))
      for (i in seq_along(plots_vol)) {
        if (plots_vol[[i]]$args$stat_name == stat_name_selected()) {
          data_vol <- plots_vol[[i]]$output[[1]]$data
        }
      }
      
      row_data <- rowData(D) %>% data.frame()
      names <- unlist(row_data[row_data$SUB_PATHWAY == name,]$name)
      data_vol$isSelected <- data_vol$name %in% names
      data_vol$text <- ifelse(data_vol$name %in% names, data_vol$name, "")
      
      t <- list(family = "sans serif", size = 14, color = toRGB("grey50"))
      
      plot <- data_vol %>%
        ggplot(aes(x = statistic, y = p.value, color = isSelected)) +
        geom_point(aes(text = name)) +
        scale_y_continuous(trans = reverselog_trans(10),
                           breaks = scales::trans_breaks("log10", function(x) 10^x),
                           labels = scales::trans_format("log10", scales::math_format(10^.x))) +
        labs(y = "p-value") +
        ggtitle(paste0(stat_name_selected(), "-", name)) +
        ggrepel::geom_text_repel(aes(label = name), max.overlaps = Inf)
      session_store$mod4_vol <- ggplotly(plot, source = "mod4_sub_vol") %>%
        layout(dragmode = "lasso",
               legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3)) %>%
        add_text(text=~data_vol$text,
                 textposition="top right",
                 showlegend = T)
      session_store$mod4_vol
    }
  })
  
  # download button
  output$mod4_download_plotly_volcano <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod4_vol), file, selfcontained = TRUE)
    }
  )
  
  # box-scatter plot
  output$mod4_box_scatter <- renderPlotly({
    d1 <- event_data("plotly_click", source = "mod4_sub_bar")
    d2 <- event_data("plotly_click", source = "mod4_sub_vol")
    if (!is.null(d1) & !is.null(d2)) {
      plots_bar <- mtm_res_get_entries(D, c("plots", "stats"))
      for (i in 1:(length(plots_bar)-1)) {
        if (plots_bar[[i]]$args$stat_list == stat_name_selected()) {
          data_bar <- plots_bar[[i]]$output[[1]]$data
        }
      }
      lvls <- rev(levels(data_bar$label))
      label <- lvls[round(as.numeric(d1$y))]
      name <- data_bar[data_bar$label == label, ]$name
      
      plots_vol <- mtm_res_get_entries(D, c("plots", "volcano"))
      for (i in seq_along(plots_vol)) {
        if (plots_vol[[i]]$args$stat_name == stat_name_selected()) {
          data_vol <- plots_vol[[i]]$output[[1]]$data
        }
      }
      
      row_data <- rowData(D) %>% data.frame()
      names <- unlist(row_data[row_data$SUB_PATHWAY == name,]$name)
      data_vol <- data_vol[data_vol$name %in% names, ]
      name2 <- data_vol[as.numeric(d2$pointNumber) + 1, ]$name[1]
      
      rd <- rowData(D) %>%
        as.data.frame() %>%
        dplyr::mutate(var = rownames(D))
      
      stat <- maplet::mtm_get_stat_by_name(D, stat_name_selected()) %>%
        dplyr::inner_join(rd, by = "var")
      
      data_box <- D %>%
        maplet:::mti_format_se_samplewise() %>%
        tidyr::gather(var, value, dplyr::one_of(rownames(D)))%>%
        dplyr::rename("Age Regression" = "Age",
                      "Comparison Outcome1" = "outcome1",
                      "Comparison Outcome2" = "outcome2",
                      "Diagnosis Analysis" = "Diagnosis") %>%
        dplyr::select(dplyr::one_of(c("var","value", stat_name_selected()))) %>%
        dplyr::inner_join(stat[,dplyr::intersect(colnames(stat),
                                                 c('var','statistic','p.value','p.adj','name'))], by = "var") %>%
        dplyr::select(-var)
      
      data_box <- data_box[data_box$name==name2, ]
      p.value <- signif(mean(data_box$p.value), 3)
      p.adj <- signif(mean(data_box$p.adj), 3)
      plot <- data_box %>%
        ggplot(aes(x = sym(stat_name_selected()),
                   y = value)) +
        geom_point(aes(text = paste0("Name:", name))) +
        geom_smooth(method = "lm", se = FALSE) +
        ggtitle(paste0(stat_name_selected(), "-", name, "-", name2))
      session_store$mod4_box <- ggplotly(plot) %>%
        layout(dragmode = "lasso") %>%
        add_annotations(
          x = min(data_box$Age) + 5,
          y = max(data_box$value),
          text = paste0("P-value: ", p.value),
          showarrow = F
        ) %>%
        add_annotations(
          x = min(data_box$Age) + 5,
          y = max(data_box$value) - 1,
          text = paste0("P.adj: ", p.adj),
          showarrow = F
        )
      session_store$mod4_box
    }
  })
  
  # download button
  output$mod4_download_plotly_scatter <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".html", sep = "")
    },
    content = function(file) {
      saveWidget(as_widget(session_store$mod4_box), file, selfcontained = TRUE)
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
                             choices = c("SUPER_PATHWAY", "SUB_PATHWAY"),
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
      ) + geom_boxplot()
      ggplotly(p)
    } else if(mod5_input()[2]==FALSE & mod5_input()[4]==TRUE) {
      p <- ggplot(as.data.frame(colData(D)),
             aes(!!sym(mod5_input()[1]), !!sym(mod5_input()[3]))
      ) + geom_boxplot()
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
      plot_ly(labels = ~var, values = ~count, type = 'pie') %>% 
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
}

# Run the application 
shinyApp(ui = ui, server = server)