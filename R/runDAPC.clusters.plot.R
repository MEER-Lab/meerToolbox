runDAPC.clusters.plot <- function(dapc_obj, gen_obj, palette = "Dark2", base_size = 15) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  
  # Extract posterior probabilities and sample IDs
  q.dat <- dapc_obj$posterior %>%
    as_tibble(rownames = "SampleID") %>%
    mutate(SiteID = adegenet::pop(gen_obj)) %>%
    pivot_longer(cols = tidyselect::matches("^[0-9]+$"),  # matches cluster columns named "1", "2", etc.
                 names_to = "cluster",
                 values_to = "prob")
  
  # Plot
  ggplot(q.dat, aes(x = SampleID, y = prob, fill = factor(cluster))) +
    geom_col(color = "gray", size = 0.1) +
    facet_grid(~fct_inorder(SiteID), switch = "x", scales = "free", space = "free") +
    labs(x = "Population", y = "Ancestry", fill = "Genetic Lineage") +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = expansion(add = 1)) +
    scale_fill_brewer(palette = palette) +
    theme_minimal(base_size = base_size) +
    theme(panel.spacing.x = unit(0.01, "lines"),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          strip.text.x = element_text(angle = -60))
}