gl.thinByDistance <- function(gl, min.dist = 10000) {
  pos <- position(gl)
  if (is.null(pos)) stop("No position data found in genlight object.")
  
  # Order loci by position
  ord <- order(pos)
  pos <- pos[ord]
  keep <- c(ord[1])
  
  for (i in ord[-1]) {
    if (pos[i] - pos[keep[length(keep)]] >= min.dist) {
      keep <- c(keep, i)
    }
  }
  
  # Subset genlight object
  gl_thinned <- gl[, keep]
  return(gl_thinned)
}
