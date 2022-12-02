merge_clusters <- function(pairs_to_merge,
                           clusters){

  for (l in 1:nrow(pairs_to_merge)) {
    pair_to_merge <- pairs_to_merge[l,]

    i         <- min(pair_to_merge)
    j         <- max(pair_to_merge)

    if(i != j){
      # merge clusters
      clusters[clusters == j] <- i
      clusters[clusters > j] <- clusters[clusters > j] - 1

      # update the rest of the table with the new clusters
      pairs_to_merge[pairs_to_merge == j] <- i
      pairs_to_merge[pairs_to_merge > j] <- pairs_to_merge[pairs_to_merge > j] - 1
    }
  }

  return(clusters)
}
