library(tidyr)
netinputsdir <- "C:/FHWA_R2/Network/CleanData"
spatialinputsdir <- "C:/FHWA_R2/spatial_boundary/CleanData"
demandinputsdir <- "C:/FHWA_R2/Demand/CleanData"
figuredir <- "C:/FHWA_R2/Network/Figures"

##### LOAD DATA #####
xwalk <- read.csv(file.path(spatialinputsdir,"cleaned_lodes8_crosswalk_with_ID.csv")) 
xwalk <- xwalk %>% dplyr::select(trct,cbsa,cbsaname,spatial_id)
colnames(xwalk)[colnames(xwalk) == "trct"] = "tract"
tracts <- read.csv(file.path(netinputsdir,"network_microtype_metrics_2.csv"))
tracts <- merge(x = tracts, y = xwalk, by = "tract") # 83,612 tracts
geotype <- read.csv(file.path(demandinputsdir,"geotype_inputs.csv")) 
geotype$est_population <- geotype$pop_density*(0.000247105*geotype$ALAND)
geotype <- geotype %>% dplyr::select(-"ALAND",-"AWATER",-"pct_water",-"laneMiles",-"total_jobs")
tracts <- merge(x = tracts, y = geotype, by = "spatial_id") # 83,612 tracts
demand <- read.csv(file.path(demandinputsdir,"microtype_inputs_demand_V2.csv")) 
colnames(demand)[colnames(demand) == "GEOID"] = "tract"
demand <- demand %>% dplyr::select(-"ALAND",)
osm <- read.csv(file.path(netinputsdir,"OSMNX/osm_metrics.csv"))
osm <- osm %>% dplyr::select(-"STATEFP",-"TRACTCE",-"X",-"ALAND")
colnames(osm)[colnames(osm) == "GEOID"] = "tract"
## Filter Non-partition tracts (>=10km^2 or not CBSA)
cbsa_df <- merge(x = tracts, y = demand, by = "tract",all.x=TRUE) # 83,612 tracts (24 without osm data, 26 without demand data); 5946 (2 without demand data) 30,653 tracts
cbsa_df <-  merge(x = cbsa_df, y = osm, by = "tract",all.x=TRUE) # 83,612 total tracts; 30,653 tracts

##### PREPARE FEATURES ######
features <- c("pct_controlp","RRS",#"avg_iri",
              "percent_laneMiles_f_system_1_2","percent_road_f_system_3","percent_road_f_system_4","percent_road_f_system_5_7",
              "lane_mile_density","lane_miles_per_capita",
              "self_loop_proportion","street_density","avg_street_length","circuity_avg","streets_per_node_avg",
              #"edge_count",#"node_count",
              "jobs_per_acre","pop_per_acre","Impervious.Developed","Developed.Open.Space","jobs_resident_bal_imp","job_sink_mag",
              "intersection_density",
              "intersection_per_street"
)
id_fields <- c("tract","spatial_id","cbsa","is_cbsa") #"tract","network_id"
colSums(is.na(cbsa_df%>% select(any_of(c(id_fields,features)))))
colSums(is.na(cbsa_df%>% select(any_of(c("populationE","pop_per_acre","land_area_acre","est_population")))))
# calculate derived columns: percent_laneMiles_f_system_1_2, lane_miles_per_capita
cbsa_df$percent_laneMiles_f_system_1_2 <- cbsa_df$percent_laneMiles_f_system_1+cbsa_df$percent_laneMiles_f_system_2
# replace 0 tract population with cbsa population estimate
cbsa_df$aland_acre <- cbsa_df$aland*0.000247105
cbsa_df$est_population <- cbsa_df$pop_density*cbsa_df$aland_acre
cbsa_df <- cbsa_df %>% mutate(populationE = replace_na(populationE, 0))
cbsa_df$populationE_imp <- cbsa_df$populationE
cbsa_df$populationE_imp[cbsa_df$populationE == 0] <- cbsa_df$est_population[cbsa_df$populationE == 0]
cbsa_df$lane_miles_per_capita <- cbsa_df$laneMiles / cbsa_df$populationE_imp

cbsa_df$est_jobs <- cbsa_df$job_density*cbsa_df$aland_acre
cbsa_df$jobs_resident_bal_imp <- cbsa_df$jobs_resident_bal
cbsa_df$jobs_resident_bal_imp[cbsa_df$populationE == 0] <- cbsa_df$est_jobs[cbsa_df$populationE == 0]/cbsa_df$populationE_imp[cbsa_df$populationE == 0]
# replace NaN job_sink_mag with 0
cbsa_df <- cbsa_df %>% mutate(job_sink_mag = if_else(is.na(job_sink_mag), 0, job_sink_mag))
colSums(is.na(cbsa_df%>% select(any_of(c(id_fields,features)))))
# impute data from aggregate by spatial id
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(RRS = ifelse(is.na(RRS), mean(RRS, na.rm = TRUE), RRS)) %>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(lane_miles_per_capita = ifelse(is.na(lane_miles_per_capita), sum(laneMiles, na.rm = TRUE) / sum(populationE_imp, na.rm = TRUE), lane_miles_per_capita))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(jobs_resident_bal_imp = ifelse(is.na(jobs_resident_bal_imp), sum(total_jobs, na.rm = TRUE) / sum(populationE_imp, na.rm = TRUE), jobs_resident_bal_imp))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(pop_per_acre = ifelse(is.na(pop_per_acre), sum(populationE_imp, na.rm = TRUE) / sum(land_area_acre, na.rm = TRUE), pop_per_acre))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(jobs_per_acre = ifelse(is.na(jobs_per_acre), sum(total_jobs, na.rm = TRUE) / sum(land_area_acre, na.rm = TRUE), jobs_per_acre))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate( streets_per_node_avg = ifelse(is.na( streets_per_node_avg), sum(streets_per_node_avg*street_count, na.rm = TRUE) / sum(street_count, na.rm = TRUE),  streets_per_node_avg))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(self_loop_proportion = ifelse(is.na(self_loop_proportion), sum(self_loop_proportion*edge_count, na.rm = TRUE)/sum(edge_count, na.rm = TRUE), self_loop_proportion))%>% ungroup()
cbsa_df$total_street_length <- cbsa_df$avg_street_length*cbsa_df$street_count
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(avg_street_length = ifelse(is.na(avg_street_length), sum(total_street_length, na.rm = TRUE)/sum(edge_count, na.rm = TRUE), avg_street_length))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(circuity_avg = ifelse(is.na(circuity_avg), sum(edge_length_total, na.rm=TRUE)/sum(edge_length_total/circuity_avg, na.rm = TRUE), circuity_avg))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(street_density = ifelse(is.na(street_density), sum(street_density*aland, na.rm = TRUE)/sum(aland, na.rm=TRUE), street_density))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(intersection_density = ifelse(is.na(intersection_density), sum(intersection_density*aland, na.rm = TRUE)/sum(aland, na.rm=TRUE), intersection_density))%>% ungroup()
cbsa_df$intersection_per_street <- cbsa_df$intersection_density/cbsa_df$street_density
cbsa_df <- cbsa_df %>% mutate(intersection_per_street = if_else(is.na(intersection_per_street), 0, intersection_per_street))
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(edge_count = ifelse(is.na(edge_count), sum(edge_count, na.rm = TRUE), edge_count))%>% ungroup()

cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(Impervious.Developed = ifelse(is.na(Impervious.Developed), sum(Impervious.Developed*aland, na.rm = TRUE)/sum(aland,na.rm=TRUE), Impervious.Developed))%>% ungroup()
cbsa_df <- cbsa_df %>%
  group_by(spatial_id) %>%
  mutate(Developed.Open.Space = ifelse(is.na(Developed.Open.Space), sum(Developed.Open.Space*aland, na.rm = TRUE)/sum(aland,na.rm=TRUE), Developed.Open.Space))%>% ungroup()

cbsa_df[sapply(cbsa_df, is.infinite)] <- 0 

colSums(is.na(cbsa_df%>% select(c(id_fields,features))))


##### RUN CLUSTERING DIAGNOSTICS #####
inputs <- cbsa_df %>% select(c(id_fields,features)) %>% distinct()
# check number of missing values for each variables, 
colSums(is.na(inputs)) # 6753 tracts have no iri values; now 2844 with no avg iri
write.csv(inputs, file.path(netinputsdir,"network_geotype_ALL_unscaled_data_byTracts.csv"), row.names=FALSE)

inputs_scaled <- cbind(as.data.frame(scale(inputs %>% dplyr::select(-is_cbsa,-spatial_id,-cbsa,-tract,) )), inputs %>% dplyr::select(tract,is_cbsa,spatial_id,cbsa))# scale to mean = 0, sd = 1, returns a matrix,

inputs_cbsa_scaled <- inputs_scaled%>% filter(is_cbsa==1)
inputs_non_cbsa_scaled <- inputs_scaled%>% filter(is_cbsa==0)
colSums(is.na(inputs_cbsa_scaled))
colSums(is.na(inputs_non_cbsa_scaled))

# general correlations between all inputs
cor_all <- cor(inputs_cbsa_scaled%>% dplyr::select(-is_cbsa,-spatial_id,-tract,-cbsa), use = "pairwise.complete.obs") #everything
colnames(cor_all)
corrplot(cor_all, tl.col = "black", tl.cex = .7, 
         method = "circle", type = "upper")

ggsave(file = file.path(figuredir,'CBSA_cluster_corr.png'),
       height = 3, width = 5) 
cor_all <- cor(inputs_non_cbsa_scaled%>% dplyr::select(-is_cbsa,-spatial_id,-tract,-cbsa), use = "pairwise.complete.obs") #everything
colnames(cor_all)
corrplot(cor_all, tl.col = "black", tl.cex = .7, 
         method = "circle", type = "upper")

ggsave(file = file.path(figuredir,'NON_CBSA_cluster_corr.png'),
       height = 3, width = 5) 

# First determine the number of clusters
maxC=10
for(which_input in ("non_cbsa","cbsa"){
  if(which_input=='cbsa'){
    this_input <- inputs_cbsa_scaled%>% dplyr::select(-is_cbsa,-spatial_id,-tract,-cbsa)
  }else{
    this_input <- inputs_non_cbsa_scaled%>% dplyr::select(-is_cbsa,-spatial_id,-tract,-cbsa)
  }
  reduced.res = NULL
  for(k in 2:maxC){
    print(paste('k=',k))
    # using reduced dimension data
    cluster.inputs <- clara(this_input, k, metric = "euclidean", #euclidean distance metric
                            stand = FALSE, samples = 1000, sampsize = 100, pamLike = TRUE) #partition around mediods
    col_name = paste0('cluster',k)
    if(which_input=='cbsa'){
      inputs_cbsa_scaled[[col_name]] <- as.factor(cluster.inputs$cluster)
    }else{
      inputs_non_cbsa_scaled[[col_name]] <- as.factor(cluster.inputs$cluster)
    }
    
    indexVal.inputs = intCriteria(as.matrix(this_input),
                                  cluster.inputs$clustering,c("Davies_Bouldin","Silhouette"))
    
    tmp = data.frame(k = k,
                     Inv_Bouldin = 1/indexVal.inputs$davies_bouldin, 
                     Silhouette = indexVal.inputs$silhouette)
    
    reduced.res = rbind(reduced.res,tmp)
    
    rm(tmp);rm(indexVal.inputs)
    
  }
  write.csv(inputs_non_cbsa_scaled, file.path(netinputsdir,"network_geotype_non_cbsa_clustered_scaled_data_exEdgeCount_byTracts_IntersectionPerStreet.csv"), row.names=FALSE)
  #write.csv(inputs_cbsa_scaled, file.path(netinputsdir,"network_geotype_cbsa_clustered_scaled_data_exEdgeCount_byTracts_IntersectionPerStreet.csv"), row.names=FALSE)
  
}



ggplot(reduced.res, aes(x = k, y = Inv_Bouldin)) +
  theme_bw()+
  geom_point(color = '#f59ff2') +
  geom_line(color = '#f59ff2') +
  geom_point(aes(x = k, y = Silhouette), color = "#7be2ed")+
  geom_line(aes(x =k, y = Silhouette),color = "#7be2ed") +
  scale_y_continuous(
    # Features of the first axis
    name = "Inverse DBI & Silhouette"
    # Add a second axis and specify its features
    #sec.axis = sec_axis(~1, name="Silhouette")
    ) + 
  theme(
    axis.title.y = element_text(color ='#f59ff2', size=13),
    axis.title.y.right = element_text(color = "#7be2ed", size=13)) +
  xlab('Number of Clusters') +
  scale_x_continuous(breaks = 2:maxC, labels = 2:maxC)+
  ggtitle('Method = PAM-LIKE')
ggsave(file = file.path(figuredir,paste(which_input,'_cluster_validity_IntersectionPerStreet.png')),
       height = 3, width = 5) 




# k=2 has highest values
this_input <- inputs_non_cbsa_scaled%>% dplyr::select(-is_cbsa,-spatial_id,-tract,-cbsa)
for(k in c(2,5,6,7)){
  print(k)
  clusters <- clara(this_input, k, metric = "euclidean", #euclidean distance metric
                   stand = FALSE, samples = 1000, sampsize = 100, pamLike = TRUE)
  col_name = paste0('cluster',k)
  inputs_non_cbsa[[col_name]] <- as.factor(clusters$cluster)
  inputs_non_cbsa_scaled[[col_name]] <- as.factor(clusters$cluster)


  table(inputs_non_cbsa[[col_name]])
}
#inputs_cbsa_scaled$network_id <- inputs_cbsa$network_id
#inputs_cbsa_scaled$is_cbsa <- inputs_cbsa$is_cbsa

write.csv(inputs_non_cbsa, file.path(netinputsdir,"network_geotype_non_cbsa_clustered_exEdgeCount_byTracts.csv"), row.names=FALSE)
write.csv(inputs_non_cbsa_scaled, file.path(netinputsdir,"network_geotype_non_cbsa_clustered_scaled_data_exEdgeCount_byTracts.csv"), row.names=FALSE)

this_input <- inputs_cbsa_scaled%>% dplyr::select(-is_cbsa,-spatial_id,-tract,-cbsa)
for(k in c(2,5,6)){
  print(k)
  clusters <- clara(this_input, k, metric = "euclidean", #euclidean distance metric
                    stand = FALSE, samples = 1000, sampsize = 100, pamLike = TRUE)
  col_name = paste0('cluster',k)
  inputs_cbsa[[col_name]] <- as.factor(clusters$cluster)
  inputs_cbsa_scaled[[col_name]] <- as.factor(clusters$cluster)
  
  
  table(inputs_cbsa[[col_name]])
}
#inputs_cbsa_scaled$network_id <- inputs_cbsa$network_id
#inputs_cbsa_scaled$is_cbsa <- inputs_cbsa$is_cbsa

write.csv(inputs_cbsa, file.path(netinputsdir,"network_geotype_cbsa_clustered_exEdgeCount.csv"), row.names=FALSE)
write.csv(inputs_cbsa_scaled, file.path(netinputsdir,"network_geotype_cbsa_clustered_scaled_data_exEdgeCount.csv"), row.names=FALSE)


