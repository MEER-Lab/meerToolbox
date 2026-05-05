runDAPC.clusters.plot <- function(
    dapc_obj,
    gen_obj,
    order_cluster = NULL,   # numeric cluster ID to order by
    plot_title = NULL,      # NEW: optional plot title
    palette = "Dark2",
    base_size = 15
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  
  # Extract posterior probabilities and sample IDs
  q.dat <- dapc_obj$posterior %>%
    as_tibble(rownames = "SampleID") %>%
    mutate(SiteID = adegenet::pop(gen_obj)) %>%
    pivot_longer(
      cols = tidyselect::matches("^[0-9]+$"),
      names_to = "cluster",
      values_to = "prob"
    )
  
  # Optional ordering by a specific DAPC cluster
  if (!is.null(order_cluster)) {
    cluster_name <- as.character(order_cluster)
    
    if (!cluster_name %in% colnames(dapc_obj$posterior)) {
      stop("order_cluster must match a numeric cluster present in dapc_obj$posterior.")
    }
    
    order_df <- dapc_obj$posterior %>%
      as_tibble(rownames = "SampleID") %>%
      select(SampleID, !!cluster_name := all_of(cluster_name)) %>%
      rename(order_value = !!cluster_name)
    
    q.dat <- q.dat %>%
      left_join(order_df, by = "SampleID") %>%
      arrange(SiteID, desc(order_value), SampleID) %>%
      mutate(SampleID = factor(SampleID, levels = unique(SampleID)))
  }
  
  # Base plot
  p <- ggplot(q.dat, aes(x = SampleID, y = prob, fill = factor(cluster))) +
    geom_col(color = "gray", linewidth = 0.1) +
    facet_grid(~fct_inorder(SiteID), switch = "x", scales = "free", space = "free") +
    labs(
      x = "Population",
      y = "Ancestry",
      fill = "Genetic Lineage"
    ) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 1)) +
    scale_fill_brewer(palette = palette) +
    theme_minimal(base_size = base_size) +
    theme(
      panel.spacing.x = unit(0.01, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid = element_blank(),
      strip.text.x = element_text(angle = -90)
    )
  
  # Add title only if provided
  if (!is.null(plot_title)) {
    p <- p + ggtitle(plot_title)
  }
  
  p
}