plot_locus_allele_depth <- function(df, locus_name, outdir) {
  
  dat <- df %>%
    filter(locus == locus_name) %>%
    pivot_longer(-locus, names_to = "sample", values_to = "counts") %>%
    separate(counts, into = c("a1", "a2"), sep = ",", convert = TRUE) %>%
    filter(!is.na(a1) & !is.na(a2))
  
  p <- ggplot(dat, aes(x = a1, y = a2)) +
    geom_point(alpha = 0.5, size = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    scale_x_continuous("Allele 1 read depth") +
    scale_y_continuous("Allele 2 read depth") +
    ggtitle(paste("Allele balance:", locus_name)) +
    theme_bw(base_size = 12)
  
  ggsave(
    filename = file.path(outdir, paste0(locus_name, ".png")),
    plot = p,
    width = 5,
    height = 5,
    dpi = 300
  )
}