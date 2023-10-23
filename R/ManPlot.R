#' This function provides a manhattan plot for the magnitude of eigenvector loadings.
#' 
#' @param V contains the eigenevectors corresponding to the loadings
#' @param predictor.names names of predictors corresponding to loadings in V
#' @param group.structure grouping structure of loadings in the Manhattan Plot
#' @param specific.groups allows certain groups to be selected
#' @param exclude.groups allows certain groups to be excluded
#' @param man.thresh draws a horizontal line at a specified threshold
#' @param special.colors data frame of group structures with corresponding color palettes
#' @param optional.index gives values if index values differ from predictor order in V
#' @param point.size size of the plot points
#' @keywords CIFTI manhattan plot
#' @export
#' @examples CifManPlot(Zlist, sepAnalysis, nresp, man.thresh, resp.names, 
#'   specific.groups=NULL, exclude.groups=NULL, special.colors=NULL)


ManPlot <- function(V,predictor.names=NULL,group.structure=NULL,specific.groups=NULL, exclude.groups=NULL,man.thresh=0, special.colors=NULL, optional.index=NULL,point.size=1.5){
  V <- as.matrix(V)
  p <- dim(V)[[1]]
  nvec <- dim(V)[[2]]
  unique.group <- unique(group.structure)
  num.group <- length(unique.group)
  
  if(is.null(predictor.names)) predictor.names <- paste0("X",c(1:p))
  
  if(is.null(optional.index)){
    index <- c(1:p)
  } else{
    if(length(optional.index)!=p) stop("Specify same number of indices as predictors")
  }
  
  if(!is.null(special.colors)){
    dfcolor <- data.frame(special.colors)
    colnames(dfcolor) <- c("group", "color")
    dfcolor$group <- as.factor(dfcolor$group)
    if(num.group!=dim(dfcolor)[[1]]) stop("Specify same number groups as colors")
  }
  
  for(i in 1:nvec){
  
  df <- data.frame(pred=predictor.names, group=group.structure, index=index)
  df$group <- as.factor(df$group)
  
  if(!is.null(specific.groups) || !is.null(exclude.groups)){
    if(!is.null(specific.groups)){
      df <- df %>% dplyr::filter(group %in% specific.groups)
    }
    if(!is.null(exclude.groups)){
      df <- df %>% dplyr::filter(!group %in% exclude.groups)
    }
  }
  
  if(!is.null(special.colors)){
    if(!is.null(specific.groups)){
      dfcolor <- dfcolor %>% dplyr::filter(group %in% specific.groups)
    }
    if(!is.null(exclude.groups)){
      dfcolor <- dfcolor %>% dplyr::filter(!group %in% exclude.groups)
    }
  }
  
    if(nvec==1){
      df$V <- abs(V)
    } else{
      df$V <- abs(V[,i])
    }
  
  df <- df %>% dplyr::group_by(group) %>%
    dplyr::arrange(index) %>%
    dplyr::mutate(newindx = 1:dplyr::n()) %>%
    dplyr::arrange(group,newindx) 
  
  df2 <-  df %>% dplyr::group_by(group) %>%
    dplyr::summarise(max_idx = max(newindx)) %>%
    dplyr::mutate(idx_add = dplyr::lag(cumsum(max_idx), default = 0)) %>% 
    dplyr::select(group, idx_add) %>%
    dplyr::ungroup() 
  
  plot_data <- df %>% dplyr::inner_join(df2, by="group") %>%
    dplyr::mutate(idx_cum = newindx + idx_add)
  
  axis_set <- plot_data %>% 
    dplyr::group_by(group) %>% 
    dplyr::summarize(center = mean(idx_cum)) %>%
    dplyr::left_join(dfcolor, by="group") %>%
    dplyr::ungroup()
  
  if(is.null(special.colors)) dfcolor <- data.frame(color=rep(c("#276FBF", "#183059"), unique(length(axis_set$group))))
  
    fname <- paste0("Manhattan_Vector_", i, ".png")
    gg_plot_title <- paste0("Manhattan Vector ", 1)
    
    # ylim <- plot_data %>% 
    # dplyr::mutate(ylim = (max(V) + 0.2*max(V))) %>% 
    # dplyr::pull(ylim)
    
    ylim <- 1.2*max(plot_data$V)
      
    if(man.thresh > 0){
      manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = V, color = as_factor(group))) +
          #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
      geom_point(alpha = 0.75, size=point.size) +
      geom_hline(yintercept=man.thresh, color="red", linetype="solid") +
      scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
      scale_color_manual(values = axis_set$color) +
      scale_size_continuous(range = c(0.25,2.5)) +
      labs(x = NULL, y = "Absolute Value \nEigenvector \nLoadings") + 
      ggtitle(gg_plot_title) +
      theme_minimal() + 
      theme( 
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = ggtext::element_markdown(),
        axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
        )
      } else{
      manhplot <- ggplot2::ggplot(data=plot_data, mapping=aes(x = idx_cum, y = V, color = as_factor(group))) +
      #geom_hline(yintercept = 2.2, color = "grey40", linetype = "dashed") + 
      geom_point(alpha = 0.75, size=point.size) +
      scale_x_continuous(label = axis_set$group, breaks = axis_set$center) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(axis_set$group)))) +
      scale_color_manual(values = axis_set$color) +
      scale_size_continuous(range = c(0.25,2.5)) +
      labs(x = NULL, y = "Absolute Value \nEigenvector \nLoadings") + 
      ggtitle(gg_plot_title) +
      theme_minimal() + 
      theme( 
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.y = ggtext::element_markdown(),
        axis.text.x = element_text(angle = 60, size = 6, vjust = 0.5)
        )
      }
    png(filename=fname)
    print(manhplot)
    dev.off()
  }
}