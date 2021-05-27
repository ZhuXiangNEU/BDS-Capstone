rm(list=ls())

library(shiny)
library(maplet)
library(tidyverse)
library(DT)
library(plotly)
library(shinyWidgets)



########################################################################## 
########################   Data Exploration  #############################
###########################################################
############### 
# load SE
load("SE.Rdata")

# define pathway annotation column (extracted from corresponding stat_bar plot)
pwvar <- "pathway"
# define threshold for significance (extracted from corresponding stat_bar plot)
alpha <- 0.2
# define stat_names (these should be extracted dynamically from the SE object
comp <- c("Age Regression",
          "Diagnosis Analysis",
          "Comparison Outcome1",
          "Comparison Outcome2")
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
                         downloadButton("download_plotly_bar", "download bar plot"),
                         plotlyOutput("mod2.bar", height = 600),
                         br(),
                         downloadButton("download_plotly_volcano", "download volcano plot"),
                         plotlyOutput("mod2.vol", height = 600),
                         br(),
                         downloadButton("download_plotly_scatter", "download scatter plot"),
                         plotlyOutput("mod2.box", height = 400),
                         br()
                     )
                 )), 
        
        
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
    
    # create reactive inputs list
    mod2_input_object <- eventReactive(input$mod2_go, 
                                       {c(input$mod2.stat,
                                          input$mod2.plot2,
                                          input$mod2.plot3)}
    )
    
    session_store <- reactiveValues()
    
    #
    output$mod2.bar <- renderPlotly({
        inputs <- mod2_input_object()
        plots <- mtm_res_get_entries(D, c("plots", "stats"))
        for (i in seq_along(plots)) {
            if (plots[[i]]$args$stat_list == inputs[1]) {
                plot <- plots[[i]]$output[[1]]
            }
        }
        session_store$bar <- ggplotly(plot, source = "sub_bar") %>%
            layout(dragmode = "lasso",
                   legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3))
        # render plotly graph
        session_store$bar
    })
    
    # 
    output$download_plotly_bar <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".html", sep = "")
        },
        content = function(file) {
            # export plotly html widget as a temp file to download.
            saveWidget(as_widget(session_store$bar), file, selfcontained = TRUE)
        }
    )
    
    output$mod2.vol <- renderPlotly({
        inputs <- mod2_input_object()
        d <- event_data("plotly_click", source = "sub_bar")
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
            session_store$vol <- ggplotly(plot, source = "sub_vol") %>%
                layout(dragmode = "lasso",
                       legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3)) %>%
                add_text(text=~data_vol$text,
                         textposition="top right",
                         showlegend = T)
            session_store$vol
        }
    })
    
    # 
    output$download_plotly_volcano <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".html", sep = "")
        },
        content = function(file) {
            saveWidget(as_widget(session_store$vol), file, selfcontained = TRUE)
        }
    )
    
    #
    output$mod2.box <- renderPlotly({
        inputs <- mod2_input_object()
        d1 <- event_data("plotly_click", source = "sub_bar")
        d2 <- event_data("plotly_click", source = "sub_vol")
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
            session_store$box <- ggplotly(plot) %>%
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
            session_store$box
        }
    })
    
    # 
    output$download_plotly_scatter <- downloadHandler(
        filename = function() {
            paste("data-", Sys.Date(), ".html", sep = "")
        },
        content = function(file) {
            saveWidget(as_widget(session_store$box), file, selfcontained = TRUE)
        }
    )
    
}


profvis::profvis(runApp(shinyApp(ui = ui, server = server)))
