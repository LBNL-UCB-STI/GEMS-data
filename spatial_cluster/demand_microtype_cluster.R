###############
## CODE FOR FIRST-STAGE CLUSTERING TO GENERATE DEMAND MICROTYPES
# IMPORT CLEANED DATASET OF RAW INPUTS
# GENERATE FACTORS WITH PCA
# CLUSTER USING CLARA
# SAVE SHAPEFILES OF MICROTYPES

# XIAODAN XU & LING JIN
# LAWRENCE BERKELEY NATIONAL LAB
# Main File: FEB 5 2024
##################################################
# Load prep files
# CODE DIRECTORY
mywd <- "C:/Users/xiaodanxu/Documents/GitHub/GEMS-data/spatial_cluster"
setwd(mywd)
source('initialization.R')
source('functions.R')

mywd <- "C:/FHWA_R2/"
setwd(mywd)

figuredir <- "Demand/Figures"
inputsdir <- "Demand/CleanData"
datadir <- "Demand/Results"
tabdir <- 'Demand/Tables'
boundarydir <- 'spatial_boundary/CleanData'



######
# LOAD DATA
##########

# demand microtype input
FHWA_data <- fread(file.path(inputsdir,"microtype_inputs_demand.csv")) %>%
  mutate(tract = as.character(GEOID), census_urban_area = as.character(census_urban_area))

# Crosswalk to CBSAs and counties
xwalk <- fread(file.path(boundarydir, 'cleaned_lodes8_crosswalk_with_ID.csv')) %>%
  mutate(tract = as.character(trct))

rownames(FHWA_data) <- FHWA_data$tract # set tract as row name for identifier
FHWA_data[sapply(FHWA_data, is.infinite)] <- NA # replace infinite with NA

na.count <- colSums(is.na(FHWA_data)) # check number of missing values for each variables
max(na.count/dim(FHWA_data)[1]) # max missing is about 0.5%


# select variable for clustering
export = as.data.frame(FHWA_data) %>%
  #filter(water == 0) %>%
  dplyr::select(c('pop_per_acre', 'jobs_per_acre', 'jobs_resident_bal', 'job_diversity', 
  'jobs 0-1.3 miles', 'jobs 1.3-3 miles', 'jobs 3-8 miles', 'jobs >8 miles', 'remote jobs', 
  'job_sink_mag', 'Impervious Developed', 'Developed Open Space', 'census_urban_area', 'tract')) # remove identifiers 


#better define column names
colnames(export) <- c('pop_per_acre', 'jobs_per_acre', 'jobs_resident_bal', 'job_diversity', 
                      'jobs_dist_bin1', 'jobs_dist_bin2', 'jobs_dist_bin3', 'jobs_dist_bin4', 'remote_jobs', 
                      'job_sink_mag', 'impervious_developed', 'developed_open_space', 'census_urban_area', 'tract')

list_of_attr <- c('pop_per_acre', 'jobs_per_acre', 'jobs_resident_bal', 'job_diversity', 
                  'jobs_dist_bin1', 'jobs_dist_bin2', 'jobs_dist_bin3', 'jobs_dist_bin4', 'remote_jobs', 
                  'job_sink_mag', 'impervious_developed', 'developed_open_space')

# standardize numeric variables, center and scale them
for (i in 1:length(colnames(export))){
  if (is.numeric(export[1, i])){
    export[, i] <- as.numeric(scale(export[, i]))
  }else{
    export[, i] <- export[, i]
  }
}
write.csv(export, file.path(datadir, "microtypes_inputs_demand_scaled.csv"), row.names = F)

data_scaled <- imputeMissings(export)

urban_data_scaled <- data_scaled %>% filter(census_urban_area == '1')
rural_data_scaled <- data_scaled %>% filter(census_urban_area == '0')

data_scaled = data_scaled %>%
  select(-tract, -census_urban_area)

urban_data_scaled = urban_data_scaled %>%
  select(-tract, -census_urban_area)

rural_data_scaled= rural_data_scaled %>%
  select(-tract, -census_urban_area)

# running results diagnose
stargazer(FHWA_data, omit.summary.stat = c("p25", "p75"),  title = "Descriptive Statistics", 
          out = file.path(tabdir,"Sum_stats_stage1_vars.tex"))

# general correlations between variables at national level 
cor_all <- cor(data_scaled, use = "pairwise.complete.obs") #everything
colnames(cor_all)
corrplot(cor_all, tl.col = "black", tl.cex = .7, order = "AOE", # order them in terms of correlations
         method = "circle", type = "upper") # label size = tl.cex

png(file = file.path(figuredir, "scree_plot_microtypes.png"))
#Parallel analysis to select the number of factors
# Parallel analysis suggests that the number of factors =  12  and the number of components =  NA 

parallel <- fa.parallel(data_scaled, fm = 'minres', # factor method = min res (does not assume normal distribution)
                        fa = 'fa',  use = 'pairwise', show.legend = TRUE) # missing data
dev.off()

# demand microtype cluster
# determine optimal urban cluster numbers
sdata = sample_n(urban_data_scaled, 15000)

png(file = file.path(figuredir, "optimal_urban_numbers_of_clusters.png"))

fviz_nbclust(sdata, cluster::clara,  method = "silhouette",  k.max = 10) +
  geom_vline(xintercept = 3, linetype = 2)
dev.off()

png(file = file.path(figuredir, "optimal_rural_numbers_of_clusters.png"))
fviz_nbclust(rural_data_scaled, cluster::clara,  method = "silhouette",  k.max = 10) +
  geom_vline(xintercept = 2, linetype = 2)
dev.off()

#urban cluster
set.seed(1) # so results are replicable each time
tic()
cluster_urban <- clara(urban_data_scaled, 3, metric = "euclidean", #euclidean distance metric
                  stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
toc()

#rural cluster
set.seed(1) # so results are replicable each time
tic()
cluster_rural <- clara(rural_data_scaled, 2, metric = "euclidean", #euclidean distance metric
                       stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
toc()

urban_data <- export %>% filter(census_urban_area == '1')
rural_data <- export %>% filter(census_urban_area == '0')

urban_data_scaled <- urban_data_scaled %>% 
  mutate(demand_microtype = as.factor(cluster_urban$clustering),
         # adding identifier back in (do this after reshape)
         tract = urban_data$tract, census_urban_area =urban_data$census_urban_area) %>%
  left_join(xwalk) # link back to county/cbsa etc location info xwalk

rural_data_scaled <- rural_data_scaled %>% 
  mutate(demand_microtype = as.factor(cluster_rural$clustering),
         # adding identifier back in (do this after reshape)
         tract = rural_data$tract, census_urban_area =rural_data$census_urban_area) %>%
  left_join(xwalk) # link back to county/cbsa etc location info xwalk

data_scaled <- rbind(urban_data_scaled, rural_data_scaled)
data_scaled <- data_scaled %>%
  mutate(census_urban_area = case_when(census_urban_area == '1' ~ 'Urban',
                         TRUE ~ 'Rural')) %>%
  mutate(demand_microtype_comb = paste(census_urban_area, demand_microtype, sep = "_")) 

write.csv(data_scaled, file = file.path(datadir,'clustering_outputs.csv'), row.names = F)


# plotting demand microtype:

# urban

spider <- c(list_of_attr, 'demand_microtype_comb')

# make dataframe of mean of each variable for each cluster (remove all factor variables)
# keep cluster6 variable as factor
radar_data <- data_scaled %>% select(spider)
radar_data <-  aggregate(radar_data[list_of_attr], list(data_scaled$demand_microtype_comb), mean)

#radar_data <- aggregate(raw_variables, list(data_scaled$cluster6), mean )
#radar_data$Group.1 <- as.numeric(radar_data$Group.1)
#colnames(radar_data) <- c("Microtype",factor_labels)

# To use fmsb package, have to add 2 lines to the dataframe: the max and min of each variable to show on the plot
colMax <- function(radar_data) sapply(radar_data, max, na.rm = TRUE) #max
colMin <- function(radar_data) sapply(radar_data, min, na.rm = TRUE) #min

# first row is the column max, second row is the column min
radar_data <- rbind(colMax(radar_data), colMin(radar_data), radar_data)
rownames(radar_data) = c("max", "min", "Rural_1", "Rural_2", "Urban_1", "Urban_2", "Urban_3")

radar_data <- radar_data  %>% select(-Group.1 )
radar_data  <- mutate_all(radar_data , function(x) as.numeric(as.character(x)))

png(file = file.path(figuredir,"spider_12fac_clara6_cluster1.png"),height = 700, width = 700)
spider1 <- radarchart(radar_data[1:3, ], axistype = 0, seg = 4, pty = 32, pdensity = NULL, pangle = 45, 
                      pcol = "#636363", pfcol="#0000B3FF", plwd=4 , cex.main = 2, #title font size
                      cglty = 3, cglwd = 1, title = "Rural 1", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, 
                      vlabels = NULL, vlcex = 1.2, caxislabels = NULL, calcex = NULL, paxislabels = NULL, palcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster2.png"),height = 700, width = 700)
spider2 <- radarchart(radar_data[c(1:2,4), ], axistype = 0, pty = 32, cex.main = 2, pcol = "#636363", pfcol="#4500FFFF" , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Rural 2", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster3.png"),height = 700, width = 700)
spider3 <- radarchart(radar_data[c(1:2,5), ], axistype = 0, pty = 32, pcol = "#636363", cex.main = 2, pfcol="#C527D8FF" , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Urban 1", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster4.png"),height = 700, width = 700)
spider4 <- radarchart(radar_data[c(1:2,6), ], axistype = 0, pty = 32, pcol = "#636363", cex.main = 2, pfcol="#FFA35CFF"  , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Urban 2", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster5.png"),height = 700, width = 700)
spider5 <- radarchart(radar_data[c(1:2,7), ], axistype = 0,  pty = 32, pcol = "#636363", cex.main = 2, pfcol="#FFF50AFF"  , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Urban 3", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

colors_urban <- c("#0000B3FF", "#4500FFFF", "#C527D8FF", "#FFA35CFF", "#FFF50AFF")
png(file = file.path(figuredir, "spider_12_feature_clara_demand.png"),height = 700, width = 700)
spider_all <- radarchart(radar_data, axistype = 0,  pty = 32, plty = 1, 
                         pcol = colors_urban, plwd=4 , palcex = 1.2,
                         cglty = 3, cglwd = 1, title = "All microtypes", 
                         maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.4)
legend(x=1, y=1, legend = rownames(radar_data[-c(1,2),]), bty = "n", pch=20 , col=colors_urban , text.col = "black", cex=1.2, pt.cex=3)
dev.off()


#plot sample CA results

california_tracts <- tracts("CA", year = 2021, cb = TRUE)
library(stringr)
data_scaled_to_plot <- data_scaled %>% 
  select(tract, demand_microtype_comb) %>% 
  mutate(GEOID = stringr::str_pad(tract, 11, pad = "0")) %>% select(-tract)

california_tracts <- california_tracts %>%
  left_join(data_scaled_to_plot, by = 'GEOID')

png(file = file.path(figuredir, "demand_cluster_ca.png"))
plot(california_tracts[,'demand_microtype_comb'],
     pal = sf.colors(5), border = NA, key.pos = 3, key.length = 0.6)
dev.off()
# plot original population density
FHWA_data <- FHWA_data %>% 
  mutate(GEOID = stringr::str_pad(GEOID, 11, pad = "0"))

FHWA_data <- FHWA_data %>% 
  rename('jobs_dist_bin1'='jobs 0-1.3 miles', 
         'jobs_dist_bin2'="jobs 1.3-3 miles", 
         'jobs_dist_bin3'= "jobs 3-8 miles" , 
         'jobs_dist_bin4'= "jobs >8 miles", 
         'impervious_developed'="Impervious Developed", 
         'developed_open_space'="Developed Open Space")

FHWA_data_to_plot <- FHWA_data %>% 
  select(GEOID, pop_per_acre, jobs_per_acre)

california_tracts <- california_tracts %>%
  left_join(FHWA_data_to_plot, by = 'GEOID')
plot(california_tracts[,'pop_per_acre'],
     pal = sf.colors, border = NA, key.pos = 3, nbreaks = 6, logz=TRUE, key.length = 0.6)

plot(california_tracts[,'jobs_per_acre'],
     pal = sf.colors, border = NA, key.pos = 3, nbreaks = 6, logz=TRUE, key.length = 0.6)
# append demand typology to raw data
FHWA_data <- FHWA_data %>%
  left_join(data_scaled_to_plot, by = 'GEOID')

write.csv(FHWA_data, file = file.path(datadir,'clustering_outputs_with_raw_data.csv'), row.names = F)

