find_clusters <- function(delta_list, p, thresh_fuse = 1e-3) {
  ndiff <- length(delta_list)
  to_merge <- logical(ndiff)
  nrms <- numeric(ndiff)
  clusters = 1:p

  for (kk in 1:ndiff) {
    nrm <- norm(delta_list[[kk]], type = "2")
    nrms[[kk]] <- nrm
    to_merge[kk] <- nrm <= thresh_fuse
  }

  pairs <- vec2tri(1:ndiff, p)

  print(cbind(pairs, nrms))

  pairs_to_merge <- pairs[to_merge, , drop = FALSE]

  if(nrow(pairs_to_merge) != 0)
    clusters <- merge_clusters(pairs_to_merge, clusters)  # merge clusters

  return(clusters)
}
