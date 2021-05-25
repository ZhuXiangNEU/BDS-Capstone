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

table_stats <- mtm_res_get_entries(D, c("stats", "univ"))
plot_stats <- mtm_res_get_entries(D, c("plots", "stats"))

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
                         y = -.2,
                         tracegroupgap = 5),
           autosize = TRUE
    )
  
}
## box_plot_reduced
mt_plots_box_scatter_reduced <-
  function (D, x, stat_name, plot_type, correct_confounder, feat_filter = p.value < 
              0.05, feat_sort = p.value, annotation = "{sprintf('P-value: %.1e', p.value)}", 
            full_info = F, text_size = 3.88, jitter = "beeswarm", 
            restrict_to_used_samples = T, ylabel = NULL, fit_line = T, 
            fit_line_se = T, ggadd = NULL, 
            num, ...) 
  {
    stopifnot("SummarizedExperiment" %in% class(D))
    if (missing(x)) 
      stop("x must be provided")
    x <- dplyr::enquo(x)
    n <- sapply(as.list(substitute(list(...)))[-1L], deparse)
    dot_args <- names(n)
    if ("metab_filter" %in% dot_args) 
      stop("You used the old MT naming convention metab_filter. Should be: feat_filter")
    if ("metab_sort" %in% dot_args) 
      stop("You used the old MT naming convention metab_sort. Should be: feat_sort")
    if ("rows" %in% dot_args) 
      stop("\"rows\" is no longer an acceptecd argument.")
    if ("cols" %in% dot_args) 
      stop("\"cols\" is no longer an accepted argument.")
    if ("manual_ylab" %in% dot_args) 
      stop("You used the old MT naming convention manual_ylab. Should be: ylabel")
    if ("manual_ylabel" %in% dot_args) 
      stop("You used the old MT naming convention manual_ylabel. Should be: ylabel")
    if ("fitline" %in% dot_args) 
      stop("You used the old MT naming convention fitline. Should be: fit_line")
    if ("fitline_se" %in% dot_args) 
      stop("You used the old MT naming convention fitline_se. Should be: fit_line_se")
    Ds <- D
    if (!missing(correct_confounder)) {
      maplet:::mti_logstatus(glue::glue("correcting for {correct_confounder}"))
      Ds <- maplet:::mti_correctConfounder(Ds, correct_confounder)
    }
    rd <- rowData(Ds) %>% as.data.frame() %>% dplyr::mutate(var = rownames(Ds))
    if (!missing(stat_name)) {
      stat <- maplet::mtm_get_stat_by_name(Ds, stat_name) %>% 
        dplyr::inner_join(rd, by = "var")
    }
    else {
      stat <- rd
      restrict_to_used_samples <- F
    }
    if (!missing(feat_filter)) {
      feat_filter_q <- dplyr::enquo(feat_filter)
      stat <- stat %>% dplyr::filter(!!feat_filter_q)
      maplet:::mti_logstatus(glue::glue("filter features: {feat_filter_q} [{nrow(stat)} remaining]"))
    }
    if (!missing(feat_sort)) {
      feat_sort_q <- dplyr::enquo(feat_sort)
      stat <- stat %>% dplyr::arrange(!!feat_sort_q) %>% dplyr::mutate(name = factor(name, 
                                                                                     levels = unique(name)))
      maplet:::mti_logstatus(glue::glue("sorted features: {feat_sort_q}"))
    }
    dummy <- Ds %>% maplet:::mti_format_se_samplewise() %>% tidyr::gather(var, 
                                                                          value, dplyr::one_of(rownames(Ds)))
    if (plot_type == "box") {
      if (restrict_to_used_samples) {
        filterto <- maplet::mtm_get_stat_by_name(Ds, stat_name, 
                                                 fullstruct = T)$samples.used
        dummy <- dummy[filterto, ]
      }
    }
    mainvar <- x %>% dplyr::quo_name()
    if (mainvar %in% colnames(dummy) == F) 
      stop(glue::glue("No column in plot data frame with name \"{mainvar}\"."))
    if (!full_info) {
      vars <- x %>% dplyr::quo_name()
      q <- dplyr::quos(...)
      if (length(q) > 0) {
        vars <- c(vars, q %>% lapply(function(x) {
          x %>% as.character() %>% gsub("~", "", 
                                        .)
        }) %>% unlist() %>% as.vector())
      }
      vars <- unique(vars)
      plottitle <- ifelse(missing(stat_name), "", stat_name)
      if (plot_type == "box") {
        mainvar <- x %>% dplyr::quo_name()
        dummy[[mainvar]] <- as.factor(dummy[[mainvar]])
        p <- dummy %>% dplyr::select(dplyr::one_of(c("var", 
                                                     "value", vars))) %>% 
          dplyr::inner_join(stat[,dplyr::intersect(colnames(stat), c("var", "statistic", "p.value", "p.adj","name"))], by = "var") %>% 
          dplyr::select(-var) %>% 
          head(num*2*300) %>%
          ggplot() + geom_boxplot(aes(x = as.factor(!!x), 
                                      y = value, ...), outlier.shape = ifelse(jitter, 
                                                                              NA, 19)) + labs(x = NULL, y = NULL) + ggtitle(plottitle)
      }
      else {
        p <- dummy %>% dplyr::select(dplyr::one_of(c("var", 
                                                     "value", vars))) %>% dplyr::inner_join(stat[, 
                                                                                                 c("var", "statistic", "p.value", 
                                                                                                   "p.adj", "name")], by = "var") %>% 
          dplyr::select(-var) %>% 
          head(num*2*300) %>%
          ggplot() + {
            if (fit_line) 
              geom_smooth(aes(x = !!x, y = value), method = "lm", 
                          se = fit_line_se, color = "black")
            else NULL
          } + geom_point(aes(x = !!x, y = value, ...)) + labs(x = dplyr::quo_name(x), 
                                                              y = "feature") + ggtitle(plottitle)
      }
    }
    else {
      if (plot_type == "box") {
        plottitle <- ifelse(missing(stat_name), "", 
                            stat_name)
        p <- dummy %>% dplyr::inner_join(stat, by = "var") %>% 
          head(num*2*300) %>%
          ggplot() + geom_boxplot(aes(x = !!x, y = value, 
                                      ...), outlier.shape = ifelse(jitter, NA, 19)) + 
          labs(x = NULL, y = NULL) + ggtitle(plottitle)
      }
      else {
        plottitle <- ifelse(missing(stat_name), "", 
                            stat_name)
        p <- dummy %>% dplyr::inner_join(stat, by = "var") %>% 
          head(num*2*300) %>%
          ggplot() + {
            if (fit_line) 
              geom_smooth(aes(x = !!x, y = value), method = "lm", 
                          se = fit_line_se, color = "black")
            else NULL
          } + geom_point(aes(x = !!x, y = value, ...)) + labs(x = dplyr::quo_name(x), 
                                                              y = "feature") + ggtitle(plottitle)
      }
    }
    if (plot_type == "box") {
      if (!is.null(ylabel)) {
        p <- p + ylab(ylabel)
      }
      else {
        r <- Ds %>% maplet::mtm_res_get_entries(c("pre", 
                                                  "trans", "log"))
        if (length(r) > 0) {
          p <- p + ylab(r[[1]]$logtxt)
        }
      }
      if (jitter == "beeswarm") {
        p <- p + ggbeeswarm::geom_beeswarm(aes(x = !!x, y = value, 
                                               ...))
      }
      else if (jitter == "jitter") {
        p <- p + geom_jitter(aes(x = !!x, y = value, ...))
      }
    }
    if (!missing(annotation)) {
      data_annotate <- stat %>% dplyr::mutate(annotate = glue::glue(annotation)) %>% 
        dplyr::distinct(name, annotate)
      p <- p + geom_text(data = data_annotate, aes(label = annotate), 
                         x = -Inf, y = Inf, hjust = -0.05, vjust = 1.05, size = text_size)
    }
    if (!is.null(ggadd)) 
      p <- p + ggadd
    if (length(unique(stat$name)) == 0) {
      p <- list(ggplot() + geom_text(aes(x = 0, y = 0, label = "no plots"), 
                                     size = 10))
      output2 <- NULL
    }
    else {
      p <- p + facet_wrap(. ~ name, scales = "free_y", 
                          ncol = 2)
      p <- list(p)
      output2 <- ceiling(length(unique(stat$name))/2)
    }
    funargs <- maplet:::mti_funargs()
    D %<>% maplet:::mti_generate_result(funargs = funargs, logtxt = sprintf("Feature ", 
                                                                            ifelse(plot_type == "box", "boxplots", "scatter plots"), 
                                                                            ", aes: %s", maplet:::mti_dots_to_str(...)), output = p, 
                                        output2 = output2)
    D
  }
## equilizer plot

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
                         plotlyOutput('mod4_stats_bar', height = 730),
                         br(), 
                         # equalizer plotly
                         plotlyOutput('mod4_equalizer', height = 730),
                         br(), 
                         # box/scatter plotly
                         plotlyOutput('mod4_box_scatter', height = 730)
               )
             )
    ),
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
  ## create dynamic choices of number of box plots to render 
  output$mod1_box_plot_num_ui <- renderUI({
    if(input$mod1_select_object=="box"){
      ## limit plots to specified stat_name
      obj_name <- subset(obj_name, V1=="plots")
      obj_name <- subset(obj_name, V2=="box")
      output_order <- obj_name %>%
        mutate(order=seq(from=1, to=n()))
      output_order <- subset(output_order, stat_name==input$mod1_select_statname)
      box_plot <- mtm_res_get_entries(D, c("plots", "box"))[[output_order$order]]
      selectInput("mod1_box_plot_num", 
                  "Number of top box plots:", 
                  choices = 1:(box_plot$output2),
                  width = "220px"
      )
    }
  })
  
  # create reactive inputs list
  mod1_input_object <- eventReactive(input$mod1_go, 
                                     {c(input$mod1_radio,
                                        input$mod1_select_statname,
                                        input$mod1_select_object,
                                        input$mod1_box_plot_num)}
  )
  
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
  # create reactive inputs list
  mod3_input_object <- eventReactive(input$mod3_go, 
                                     {c(input$mod3_select_plot, 
                                        input$mod3_select_colData,
                                        input$mod3_scale_data,
                                        input$mod3_checkbox_factor,
                                        input$mod3_pca_data_type,
                                        input$mod3_umap_n_neighbors)}
  )
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
        (as.numeric(mod1_input_object()[4]))*150
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
        if(plots[[1]]$fun[2]=="box"&plots[[1]]$fun[3]=="scatter"&!is.null(plots[[row_n]]$output2)){
          box_plots <- D %>% 
            mt_plots_box_scatter_reduced(stat_name = mod1_input_object()[2],
                                         x = mod1_input_object()[2],
                                         plot_type = "scatter",
                                         feat_filter = p.adj < 0.1, 
                                         feat_sort = p.value,
                                         annotation = "{sprintf('P-value: %.2e', p.value)}\nP.adj: {sprintf('%.2e', p.adj)}",
                                         num=as.numeric(mod1_input_object()[4])) %>%
            mtm_res_get_entries(c("plots", "box"))
          
          box_plots[[output_order$order]]$output
        } else {
          plots[[row_n]][col_n]
        }
        
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
    table <- mtm_res_get_entries(D, c(mod1_input_object()[1], mod1_input_object()[3]))[[output_order$order]]$output$table %>%
      ## scientific notation
      mutate(statistic=formatC(statistic, format = "E", digits = 2),
             p.value=formatC(p.value, format = "E", digits = 2),
             p.adj=formatC(p.adj, format = "E", digits = 2)
      )
    
    ## put interested columns ahead
    table <- if ('term' %in% names(table)) {
      table %>%
        select(var, statistic, p.value, p.adj, term, everything())
    } else {
      table %>%
        select(var, statistic, p.value, p.adj, everything())
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
  
  
  # create reactive inputs list
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
  v <- reactiveValues()
  v$s <- NULL
  observe({
    if(!is.null(input$mod4_table_rows_selected)){
      v$s <- input$mod4_table_rows_selected
    }
  })
  
  # render statsBar plot in Mod4
  output$mod4_stats_bar <- renderPlotly({
    for (i in 2:length(plot_stats)) {
      stat_name_selected <- mod4_metabolite_table() %>%
        slice(v$s) %>%
        select(`stat name`)
      
      if (plot_stats[[i-1]]$args$stat_list == stat_name_selected) {
        plot <- plot_stats[[i-1]]$output[[1]]
      }
    }
    ggplotly(plot, source = "sub_bar") %>% 
      layout(dragmode = "lasso")
  })
  # render equilizer plot in Mod4
  output$mod4_equalizer <- renderPlotly({
    d <- event_data("plotly_click", source = "sub_bar")
    if (!is.null(d)) {
      ## select a stat_name by clicking one row of stat table
      stat_name_selected <- mod4_metabolite_table() %>%
        slice(v$s) %>%
        select(`stat name`) %>%
        as.character
      ## select a pathway by clicking one bar of statsBar plot
      sub_pathway_selected <- plot_stats[[i-1]]$output[[1]]$data %>%
        filter(label==d$x) %>%
        select(name) %>%
        as.character
      biochemical_selected <- rowData(D) %>% as.data.frame %>%
        filter(SUB_PATHWAY==sub_pathway_selected) %>%
        select(BIOCHEMICAL) %>%
        as.character
      pathway_selected <- rd %>%
        filter(BIOCHEMICAL==biochemical_selected) %>%
        select(pathway_name) %>%
        as.character
      ## equilizer plot
      eqplot$stat_name_selected$pathway_selected
    }
  })
}

# Run the application 
shinyApp(ui = ui, server = server)