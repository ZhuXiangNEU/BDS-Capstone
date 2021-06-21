rm(list=ls())

library(shiny)
library(maplet)
library(tidyverse)
library(DT)
library(plotly)
library(shinyWidgets)
library(RColorBrewer)



########################################################################## 
########################   Data Exploration  #############################
########################################################################## 
# load SE
load("SE.Rdata")

# assign an object of all stats tables
table_stats <- mtm_res_get_entries(D, c("stats", "univ"))
# assign an object of all stats plots
plot_stats <- mtm_res_get_entries(D, c("plots", "stats"))

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
reverselog_trans <- function(base = exp(1)){
    trans <- function(x) -log(x, base)
    inv <- function(x) base^(-x)
    scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                      scales::log_breaks(base = base),
                      domain = c(1e-100, Inf))
}


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
                         style = "overflow-y: auto; position: absolute; left: 25%",
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
        
        
        # Module 3
        tabPanel("Module 3"),
        
        
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
        
        
        # Module 5
        tabPanel("Module 5"),
        
        
        # Module 6
        tabPanel("Module 6")
        
        
    )
)


server <- function(input, output, session) {
    
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
            ggplot(aes(x = statistic, y = p.value, color = !!sym(legend_name), label = name)) +
            geom_point() +
            scale_y_continuous(trans = reverselog_trans(10),
                               breaks = scales::trans_breaks("log10", function(x) 10^x),
                               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
            labs(y = "p-value") +
            ggtitle(paste0(inputs[1])) +
            scale_color_brewer(palette="Dark2")
        
        if (inputs[2] == "bar") {
            plot <- plot +
                geom_point(size = 3)
                ggtitle(paste0(sub_pathway_name, "-", inputs[1]))
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
                # set the column curveNumber by color legend
                p.adj.significant <- alpha
                data_vol[, "curveNumber"] <- ifelse(data_vol$p.adj < as.numeric(p.adj.significant), 1, 0)
                data_vol_true <- data_vol[data_vol$curveNumber==1, ]
                data_vol_false <- data_vol[data_vol$curveNumber==0, ]
                # By using click info (curveNumber & ponitNumber) to get the metabolite name
                metabolite <- ifelse(d.vol$curveNumber == 1,
                                     data_vol_true[d.vol$pointNumber + 1, ]$var[1],
                                     data_vol_false[d.vol$pointNumber + 1, ]$var[1])
                term <- data_vol$term[1]
                
            }
        }
        
        # Filter the data by metabolite name
        data <- data[data$var == metabolite, ]
        
        # Treat as categorical or not?
        if (input$mod2.categorical) {
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
                ggplot(aes(x = var, y = !!sym(term))) +
                geom_boxplot() +
                #geom_jitter(width = 0.2) +
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
                              filter(var == input$mod4_metabolite) %>%
                              rename("name" = var)
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
        highlight_point <- data_vol[data_vol$isSelected==TRUE, ]
        
        plot <- data_vol %>%
            ggplot(aes(x = statistic, y = p.value, color = isSelected, label = name)) +
            geom_point() +
            geom_point(data=highlight_point, size = 3) +
            scale_y_continuous(trans = reverselog_trans(10),
                               breaks = scales::trans_breaks("log10", function(x) 10^x),
                               labels = scales::trans_format("log10", scales::math_format(10^.x))) +
            labs(y = "p-value") +
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
                    ggplot(aes(x = !!sym(term), y = var)) +
                    geom_point() +
                    ggtitle(metabolite)
            } else {
                plot <- data %>%
                    ggplot(aes(x = var, y = !!sym(term))) +
                    geom_boxplot() +
                    #geom_jitter(width = 0.2) +
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
    
}


shinyApp(ui = ui, server = server)
