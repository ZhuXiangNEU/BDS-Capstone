
# Extract the object names from result SE "D"------------
get_obj_name <- function(D){
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
  obj_name$stat_name <- 
    ifelse(is.na(obj_name$stat_name), 
           "(no stat_name)", 
           obj_name$stat_name)
  
  # count number of objects in the list of output of each metadata
  count_list <- c()
  for (i in seq_along(metadata(D)$results)){
    count_list[i] <- length(metadata(D)$results[[i]]$output)
  }
  obj_name$cnt <- count_list
  return(obj_name)
}


# Get the data set -----------------------
get_data_by_name <- function(D, args.name, plot.nmae, stat.name) {
  plots <- mtm_res_get_entries(D, c("plots", plot.nmae))
  for (i in seq_along(plots)) {
    if (plots[[i]]$args[args.name] == stat.name) {
      data <- plots[[i]]$output[[1]]$data
    }
  }
  data
}


# create reverselog_trans for log10-SCALE in volcano plot-----------------
reverselog_trans <- function (base = exp(1)){
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  scales::trans_new(paste0("reverselog-", format(base)), trans, inv,
                    scales::log_breaks(base = base),
                    domain = c(1e-100, Inf))
}


# define PCA output function for mod3 referring 'mt_plots_pca'------------
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


# define UMAP output function for mod3 referring 'mt_plots_umap'---------------
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
              do.call(aes_string, as.list(structure(c("x","y",color,hover), names = c("x","y","color","hover"))))
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


# define boxplot function in module 5---------------
mod5_boxplot <- function(D, x, x_cate, y, y_cate, fill, hover, ...){
  df <- data.frame(colData(D))
  ## categorize variable if user think it's not continuous
  if(x_cate==FALSE){
    df[, x] <- factor(df[, x])
  }
  if(y_cate==FALSE){
    df[, y] <- factor(df[, y])
  }
  # to reproduce jitter plot
  set.seed(4017)
  p <- ggplot(df,
              do.call(aes_string, as.list(structure(c(x, y, fill, hover), names = c("x","y", "fill", "hover"))))
  ) + 
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = .2, alpha = 0.5) +
    theme(legend.title = element_blank()) +
    ggtitle("Boxplot with ignored outliers and Jitter with Set Seed")
  
  fig <- ggplotly(p, tooltip=c("x", "y", "hover")) %>%
    layout(legend = list(orientation = 'h', xanchor = "center", x = 0.5, y = -0.3))
  fig
}


# define scatterplot function in mod5------------------
mod5_scatter <- function(D, x, y, hover, ...){
  df <- data.frame(colData(D))
  p <- ggplot(df,
              do.call(aes_string, as.list(structure(c(x, y, hover), names = c("x","y", "hover"))))
  ) + geom_point()
  
  ggplotly(p)
}


# define barplot function in mod5--------------------
mod5_barplot <- function(D, x, fill, hover, ...){
  df <- data.frame(colData(D))
  p <- ggplot(df,
              do.call(aes_string, as.list(structure(c(x, fill, hover), names = c("x","fill", "hover"))))
  ) + geom_bar()
  
  ggplotly(p) %>%
    layout(legend = list(orientation = "h",   # show entries horizontally
                         xanchor = "center",  # use center of legend as anchor
                         x = 0.5, ## set position of legend
                         y = -0.2,
                         tracegroupgap = 5),
           autosize = TRUE
    )
}

# SE generating function for Mod6
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