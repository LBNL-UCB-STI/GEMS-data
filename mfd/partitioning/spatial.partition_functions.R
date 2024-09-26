# find the optimum cluster
find.opt.k <- function(res){
  # function to find the optimal number
  # input res is the dataframe[k, dbi], k start from 2
  difdbi = diff(res$DBI)
  opt.k = length(difdbi) # the default is the max
  for(i in 1:(length(difdbi)-1)){
    # if the decent is greater than the rest and DBI value is among the 3 smallest
    if(res$DBI[i+1] %in% (sort(res$DBI))[1:2] & difdbi[i]<0 & (-difdbi[i]> max(-difdbi[(i+1):length(difdbi)]))){
      return(i+2)
    }
  }
  return(opt.k+2)
  
}

# compute validity metrics for a range of k value
testK <- function(cluster.df, k.range = 2:min(20,ceiling(nrow(cluster.df)/3)), method = 'kmeans'){
  # method can be kmeans or pam
  # returns the inverse dbi and asw index as a function of k
  require(clusterCrit) # for cluster internal and external validity metrics
  require(cluster)
  res = NULL
  for(k in k.range){
    print(paste('k=',k))
 #   set.seed(10) # so results are replicable each time
    if(method == 'kmeans'){
      set.seed(12345678)
      reduced <-kmeans(cluster.df, k, iter.max = 1000,nstart=100) 
      res=cbind(res,reduced$cluster)
      
    }
    
    if(method == 'pam'){
      set.seed(12345678)
      reduced = pam(cluster.df,k,stand = F, metric = "euclidean", nstart=100) # do not standardize, if need to standardize, needs to be applied in the input dataframe before calling the function
      reduced$cluster = reduced$clustering
      res=cbind(res,reduced$cluster)
      
    }
  }
    colnames(res)<-paste0('k',k.range)
    return(res)
}
 

# function to perform initial partitioning of a given continous region
# note that the region is large enough to warrant partitioning

init_partition<- function(regid){
  require(clusterCrit)
  require(rgeoda)
  require(dplyr)
  
  cdata <- tracts1 %>% ## pipe operator %>% is used to express a sequence of multiple operations in an elegant way.
    filter(regID == regid) ## filter() function is used to subset a data frame, retaining all rows that satisfy your conditions. 
  
  qwts = queen_weights(cdata) # needs to have all tracts connecting to each other
  cdata1 = cdata[paste0('hh',7:19)] # only get the columns into the clustering
  # compute validation metrics
  maxK = min(15, max(ceiling(dim(cdata1)[1]/3),4)) # need to have at least k = 2,3,4 in order for the opt.k function to work 
  #maxK = 10 
  res = data.frame(k = 2:maxK, DBI=NA)
  for(kkk in 2:maxK){
    cdat_clusters <- skater(kkk, qwts, cdata1, scale_method = 'raw', distance_method = 'euclidean') # it is important to specify raw data without normalization
    cdata$cluster = cdat_clusters$Clusters
    res$DBI[res$k == kkk] = 1/as.numeric(intCriteria(as.matrix(st_drop_geometry(cdata1)),
                                                     cdata$cluster,c("Davies_Bouldin"))) ##calculate DBI
    print(length(unique(cdata$cluster)))
    
  }
  
  ##find.opt.k(res) ## I commented this..
  # set.seed(395827)
  opt.k <- find.opt.k(res)
  cdat_clusters <- skater(opt.k, qwts, cdata1, scale_method = 'raw', distance_method = 'euclidean') # it is important to specify raw data without normalization

  cdata$cluster = cdat_clusters$Clusters 
  
  # rename the cluster in ascending order according to the density at average most congested hour
  mm = cdata[c(paste0('hh',maxHH),'cluster')] %>% st_drop_geometry() 
  names(mm)<-c('density','cluster')
  mm = mm%>%
    group_by(cluster)%>%
    summarise(density = mean(density,na.rm = T))%>%
    mutate(cluster.reorder = dense_rank(desc(-density))) %>%
    dplyr::select(-density)
  
  cdata = cdata %>%
    left_join(mm) %>%
    dplyr::select(-cluster)%>%
    rename(cluster.tmp = cluster.reorder) %>%
    mutate(cluster.tmp = str_pad(cluster.tmp, 2, pad = "0")) %>% ## padding a string
    dplyr::select(GEOID, cluster.tmp)%>%
    st_drop_geometry()
  
  tracts1 <<- tracts1 %>%
    left_join(cdata) %>%
    mutate(cluster = case_when(!is.na(cluster.tmp) ~ cluster.tmp,
                               TRUE ~ cluster)) %>%
    dplyr::select(-cluster.tmp)
  
  
  # reg <- get_tiles(tracts1, provider = "Esri.WorldTopoMap", zoom = 13)
  # 
  # # plot the microtype coverage
  # ggplot() +
  #   geom_spatraster_rgb(data = reg) +
  #   geom_sf(tracts1 %>% filter(regID == regid), mapping=aes(geometry=geometry, fill= as.character(cluster)), color='grey', size=0.1, alpha = 0.5) +
  #   coord_sf(crs = 3857)+
  #   #  scale_fill_viridis_c(option='A') +
  #   #  scale_fill_fermenter(n.breaks = length(unique(cdata$cluster)), palette = "Oranges" )+
  #   theme_bw()+
  #   ggtitle(city.sel)
  
}


## split the candidate clusters into 2 partitions
split2_partition<- function(idx){
  require(clusterCrit)
  require(rgeoda)
  require(dplyr)
  
  cdata <- tracts1 %>%
    filter(unitedID == idx) 
  
  qwts = queen_weights(cdata) ## Queen contiguity weights # needs to have all tracts connecting to each other
  cdata1 = cdata[paste0('hh',7:19)] # only get the columns into the clustering
 
  set.seed(29405)
  cdat_clusters <- skater(2, qwts, cdata1, scale_method = 'raw', distance_method = 'euclidean') # it is important to specify raw data without normalization
  
  #cdat_clusters <- skater(50, qwts, cdata1,scale_method = 'raw', distance_method = 'euclidean') # it is important to specify raw data without normalization
  cdata$cluster.tmp = str_pad(cdat_clusters$Clusters, 2, pad = '0')
  cdata = cdata %>%
    dplyr::select(GEOID, cluster.tmp)%>%
    st_drop_geometry()%>%
    distinct()

  tracts1 <<- tracts1 %>%
    left_join(cdata) %>%
    mutate(sub.cluster = case_when(
      !is.na(cluster.tmp) ~ paste(sub.cluster,cluster.tmp, sep='-'),
                               TRUE ~ sub.cluster)) %>%
    dplyr::select(-cluster.tmp)
}

  
# 
#   reg <- get_tiles(tracts1, provider = "Stamen.TonerLite", zoom = 13)
# 
#   # plot the microtype coverage
#   ggplot() +
#     geom_spatraster_rgb(data = reg) +
#     geom_sf(tracts1 %>% filter(regID == regid), mapping=aes(geometry=geometry, fill= as.character(cluster)), color='grey', size=0.1, alpha = 0.5) +
#     coord_sf(crs = 3857)+
#     #  scale_fill_viridis_c(option='A') +
#     #  scale_fill_fermenter(n.breaks = length(unique(cdata$cluster)), palette = "Oranges" )+
#     theme_bw()+
#     ggtitle(city.sel)


#-------------------------------------------------------------------------------

compute_cluster_means <- function(df){
  # Ensure 'cluster' is a factor for aggregation purposes
  df$cluster <- as.factor(df$cluster)
     
  # Define the columns of interest
  cols_of_interest <- c(paste0("hh", 7:19), "cluster")
     
  # Filter the dataframe to include only the columns of interest
  df_filtered <- df[, cols_of_interest, drop = FALSE]
  df_filtered <- df_filtered %>% st_drop_geometry()
     
  # Compute the mean feature vectors for each cluster
  mean_features <- aggregate(. ~ cluster, data = df_filtered, FUN = mean, na.rm = TRUE)
     
  # Compute the number of observations in each cluster directly
  observation_counts <- as.data.frame(table(df_filtered$cluster))
  names(observation_counts) <- c("cluster", "num_observations")
     
  # Merge mean features with observation counts
  result_df <- merge(mean_features, observation_counts, by = "cluster")
     
  return(result_df)
}
   

mini_batch_kmeans <- function(observations_df, centroids_df) {
  # Function to compute Euclidean distance
  compute_distance <- function(observation, centroid) {
    sqrt(sum((observation - centroid) ^ 2))
  }
  
  # Assign observations to the nearest cluster
  assign_clusters <- function(observation, centroids) {
    distances <- apply(centroids[, -c(1, ncol(centroids))], 1, function(centroid) compute_distance(observation, centroid))
    return(which.min(distances))
  }
  
  # Update centroids based on new assignments
  update_centroids <- function(cluster_assignments, observations, centroids) {
    
    feature_cols <- names(observations)[which(names(observations) %in% paste0("hh", 7:19))]
    
    for (cluster_id in unique(cluster_assignments)) {
      
      cluster_indices <- which(cluster_assignments == cluster_id)
      cluster_observations <- observations[cluster_indices, feature_cols, drop = FALSE]
      
      if (nrow(cluster_observations) > 0) {
        cluster_id_s = sprintf("%02d", cluster_id)
        # Convert to numeric
        cluster_observations <- data.frame(lapply(cluster_observations, as.numeric))
        
        # Calculate weighted sum of existing centroid and new observations
        existing_centroid <- centroids[centroids$cluster == cluster_id_s, feature_cols]
        existing_count <- centroids[centroids$cluster == cluster_id_s, "num_observations"]
        new_sum <- colSums(cluster_observations, na.rm = TRUE)
        total_count <- existing_count + nrow(cluster_observations)
        
        # New centroid calculation
        new_centroid <- (existing_centroid * existing_count + new_sum) / total_count
        centroids[centroids$cluster == cluster_id_s, feature_cols] <- new_centroid
        centroids[centroids$cluster == cluster_id_s, "num_observations"] <- total_count
      }
    }
    
    return(centroids)
  }
  
  # Extract feature columns from observations
  observation_features <- observations_df[, paste0("hh", 7:19), drop = FALSE]
  # Ensure observations are numeric
  observation_features <- data.frame(lapply(observation_features, as.numeric))
  
  # Assign each observation to a cluster
  cluster_assignments <- apply(observation_features, 1, function(obs) assign_clusters(obs, centroids_df))
  
  # Update centroids with new observations
  updated_centroids <- update_centroids(cluster_assignments, observations_df, centroids_df)
  
  # Create dataframe with 'GEOID' and assigned cluster
  clustered_observations <- data.frame(GEOID = observations_df$GEOID, cluster = cluster_assignments)
  
  return(list(clustered_observations = clustered_observations, updated_centroids = updated_centroids))
}
   
   

replace_with_cluster_means <- function(df) {
  # Check if df is an sf object and convert it to a standard dataframe if necessary
  if ("sf" %in% class(df)) {
    df <- as.data.frame(df)
  }
  
  # Select columns from 'hh7' to 'hh19' and group by 'cluster'
  cluster_means <- df %>%
    select(cluster, hh7:hh19) %>%
    group_by(cluster) %>%
    summarise(across(hh7:hh19, \(x) mean(x, na.rm = TRUE)), .groups = 'drop') # Ensure group structure is dropped after summarise
  
  # Join the cluster means back to the original dataframe
  df_updated <- df %>%
    select(-matches("^hh(1[0-9]|7|8|9)$")) %>%
    left_join(cluster_means, by = "cluster")
  
  return(df_updated)
}


format_clusters <- function(dataFrame) {
  dataFrame %>%
    mutate(cluster = sprintf("%02d", cluster)) # Format cluster numbers
}


## Calculate L2 distance between an orphan and neighboring cluster
L2_distance <- function(id1, id2, data){
  
  data <- st_drop_geometry(data)
  data[paste0("hh", 7:19)] <- lapply(data[paste0("hh", 7:19)], as.numeric)
  vec1 <- unique(data.matrix(subset(data, GEOID %in% id1)[paste0("hh", 7:19)])[,-c(14)])
  cluster.id <- subset(data, GEOID %in% id2)$sub.cluster
  all.vectors <- data.matrix(subset(data, sub.cluster %in% cluster.id)[paste0("hh", 7:19)])[,-c(14)]
  if (NCOL(all.vectors) > 1){ 
    vec2 <- colMeans(all.vectors)
  } else {
    vec2 <- all.vectors
  }
  dist <- as.numeric(dist(rbind(vec1, vec2)), method="euclidean")
  return(dist)
}


## Calculate DBI of the final partitioning
Davies_Bouldin_1 <- function(data){
  data <- data[c(paste0("hh", 7:19), 'unitedID')]
  code <- unique(data$unitedID)
  new_code <- c()
  m <- 1
  for (i in code){
    new_code[i] <- m
    m <- m + 1
  }
  for (j in 1:length(data$unitedID)){
    data$unitedID[j] <- as.numeric(new_code[data$unitedID[j]])
  }
  uid <- as.numeric(data$unitedID)
  data <- subset(data, select = -c(unitedID))
  DBI <- mean(index.DB(as.matrix(st_drop_geometry(data)), uid, d=NULL, centrotypes="centroids", p=2, q=2)$S)
  return(DBI)
}


## Calculate DBI
Davies_Bouldin_2 <- function(data){
  data <- data[c(paste0("hh", 7:19), 'unitedID.new')]
  code <- unique(data$unitedID.new)
  new_code <- c()
  m <- 1
  for (i in code){
    new_code[i] <- m
    m <- m + 1
  }
  for (j in 1:length(data$unitedID.new)){
    data$unitedID.new[j] <- as.numeric(new_code[data$unitedID.new[j]])
  }
  uid <- as.numeric(data$unitedID.new)
  data <- subset(data, select = -c(unitedID.new))
  DBI <- mean(index.DB(as.matrix(st_drop_geometry(data)), uid, d=NULL, centrotypes="centroids", p=2, q=2)$S)
  return(DBI)
}


# Find candidate partitions for merging orphan tracts
find_candidates <- function(dft, area.thd){
  
  s.partitions.pool <- subset(dft, cl.Ntracts %in% 1)
  s.partitions <- s.partitions.pool[s.partitions.pool$cl.aland_km2 <= area.thd,]
  
  ## Get the row index of orphan partitions
  idx <- which(dft$GEOID %in% s.partitions$GEOID)
  
  #find neighbors
  dft2 <- st_as_sf(dft)
  nb <- poly2nb(dft2, queen = TRUE) 
  neighbors <- lapply(nb, function(x) dft2$GEOID[x])
  
  ## Identify neighboring clusters for each orphan tract and store L2 distance between them
  sp.ids <- c()
  cand.ids <- c()
  members <- c()
  for (i in idx){
    if (length(neighbors[[i]]) > 0){
      sp.ids <- append(sp.ids, as.character(dft$GEOID[i]))
      cand.ids <- append(cand.ids, as.character(neighbors[[i]]))
      members[[as.character(dft$GEOID[i])]] <- as.character(neighbors[[i]])
    }
  }
  
  ## Get the set of candidates
  cand.ids <- as.vector(unique(cand.ids))
  cand.ids <- cand.ids[!cand.ids %in% s.partitions.pool$GEOID] #optional: uncomment to remove orphans from the list of candidates
  
  return(list(sp.ids, cand.ids, members))
}



# Define the new_distance function in R
new_distance <- function(ids, neighbor_ids, geoid, data) {
  # Initialize the vector with Inf
  vec <- rep(Inf, length(ids))
  
  # Vectorized operation to find indices of neighbor_ids in ids
  neighbor_indices <- which(ids %in% neighbor_ids)
  
  # Calculate L2_distance for neighbor_ids
  vec[neighbor_indices] <- sapply(ids[neighbor_indices], function(i) L2_distance(i, geoid, data))
  
  return( unlist(vec))
}




# Function to remove a row and ensure the result is always a matrix
remove_row <- function(matrix, row) {
  if (nrow(matrix) > 2) {
    new_matrix <- matrix[-row, ]
  } else {
    new_matrix <- matrix(matrix[-row, ], nrow = 1)
  }
  return(new_matrix)
}



## Iteratively merge orphan tracts
merge_orphans <- function(tracts2, area.thd){
  
  x <- 10000 ## Set arbitrary high value
  dft <- tracts2
  n <- 0
  t <- 0
  
  while (TRUE){ 
    
    if (n == 0){
      output0 <- find_candidates(dft, area.thd)
      sp.ids <- output0[1][[1]]
      cand.ids <- output0[2][[1]]
      cand.ids <- c(na.omit(cand.ids))
      members <- output0[3][[1]]
      
      if (length(cand.ids) > 0 & length(sp.ids) > 0){
        ## Create a distance matrix to store L2 distance between orphans and candidates
        dist.mat0 <- matrix(Inf, nrow = length(sp.ids), ncol = length(cand.ids))
        
        for (i in 1:length(sp.ids)){
          for (j in 1:length(cand.ids)){
            if (cand.ids[j] %in% members[[sp.ids[i]]] == TRUE){
              dist.mat0[i,j] <- L2_distance(sp.ids[i],cand.ids[j], dft)}
          }
        }
        
        dist.vec0 <- c(dist.mat0)
        dist.vec0 <- sort(unique(dist.vec0))
        dist.vec0 <- dist.vec0[dist.vec0 != c(0)]
        dist.vec0 <- dist.vec0[dist.vec0 != c(Inf)]
        dist.mat <- dist.mat0
      }
    }
    
    if (length(cand.ids) == 0 | length(sp.ids) == 0){break}
    
    #print(dim(dist.mat))
    
    #print(cand.ids)
    
    if (n > 0){ 
      sp.idx <- which(sp.ids == geoid.1)
      cand.idx <- which(cand.ids == geoid.2)
      #dist.mat <- dist.mat[-sp.idx, ]
      dist.mat <- remove_row(dist.mat, sp.idx)
      sp.ids <- sp.ids[sp.ids != geoid.1]
      neighbor.ids1 <- members[geoid.1]
      neighbor.ids2 <- members[geoid.2]
      
      vec1 <- new_distance(sp.ids, neighbor.ids1, geoid.1)
      vec2 <- new_distance(sp.ids, neighbor.ids2, geoid.2)
      
      dist.mat[,cand.idx] <- vec2
      dist.mat <- cbind(dist.mat, vec1) 
      cand.ids <- c(cand.ids, geoid.1)
    }  
    
    #match <- c()
    dist.vec <- c(dist.mat)
    dist.vec <- sort(unique(dist.vec))
    dist.vec <-  dist.vec[dist.vec != c(0)]
    dist.vec <-  dist.vec[dist.vec != c(Inf)]
    
    if (length(dist.vec) == 0){break} ## General breaking point for full convergence
    
    for (v in dist.vec){
      #print(t)
      if (t == length(dist.vec0)){break} ## Guaranteed early breaking point (after one cycle of matching sp.ids)
      t <- t + 1
      pos <- which(dist.mat == v, arr.ind = TRUE)
      for (d in 1:dim(pos)[1]){
        geoid.1 <- sp.ids[pos[d,1]] 
        geoid.2 <- cand.ids[pos[d,2]] 
        area <- subset(dft, GEOID %in% geoid.1)$cl.aland_km2 + subset(dft, GEOID %in% geoid.2)$cl.aland_km2
        n.tracts <- subset(dft, GEOID %in% geoid.2)$cl.Ntracts + 1
        #print(c(area, n.tracts))
        if (area < 4 | n.tracts < 3){break}
      }
      if (area < 4 | n.tracts < 3){break}
    }
    
    #if (area >= 4 & n.tracts >= 3){break} ## Early breaking point (after one cycle of matching sp.ids)
    if (t == length(dist.vec0)){break} 
    
    #match[as.character(geoid.1)] <- geoid.2
    s.row <- dplyr::filter(dft, GEOID %in% geoid.1)
    c.row <- dplyr::filter(dft, GEOID %in% geoid.2)
    
    c.unitedID.new <- c.row$unitedID.new
    c.cluster <- c.row$cluster
    c.sub.cluster <- c.row$sub.cluster
    
    new.cl.aland_km2 <- c.row$cl.aland_km2 + s.row$cl.aland_km2
    new.cl.Ntracts <- c.row$cl.Ntracts + s.row$cl.Ntracts
    new.density <- colMeans(st_drop_geometry(rbind(c.row[c(paste0("hh", 7:19))], s.row[c(paste0("hh", 7:19))])))
    
    member.id = c(dplyr::filter(dft, unitedID.new %in% c.unitedID.new)$GEOID, geoid.1)
    
    ## Update tracts2
    for (id in member.id){
      tracts2$unitedID.new[tracts2$GEOID == id] <- c.unitedID.new
      tracts2$cluster[tracts2$GEOID == id] <- c.cluster
      tracts2$sub.cluster[tracts2$GEOID == id] <- c.sub.cluster
      tracts2$cl.aland_km2[tracts2$GEOID == id] <- new.cl.aland_km2
      tracts2$cl.Ntracts[tracts2$GEOID == id] <- new.cl.Ntracts
      tracts2[tracts2$GEOID == id, c(paste0("hh", 7:19))] <- as.list(new.density)
    }
    
    # output.new <- find_candidates(dft, area.thd)
    # cand.ids.new <- output.new[2][[1]]
    # if (identical(cand.ids, cand.ids.new) || n >= 1000){
    #   break
    # }
    
    dft <- tracts2
    n <- n + 1
  }
  
  return(list(n, tracts2, dft))
}


   
   
   
   
   
   
   
   
   