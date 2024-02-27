###############
## CODE FOR FIRST-STAGE CLUSTERING TO GENERATE MICROTYPES
        # IMPORT CLEANED DATASET OF RAW INPUTS
        # GENERATE FACTORS WITH EXPLORATORY FACTOR ANALYSIS
        # CLUSTER USING CLARA
        # SAVE SHAPEFILES OF MICROTYPES

# NATALIE POPOVICH & LING JIN
# LAWRENCE BERKELEY NATIONAL LAB
# Main File: SEP 2 2020
# LJ 11/11/2020: use alternative set of inputs for stage 1
# 11/18/2020: decided to use 12 factors
# 11/19/2020: remove tracts with missing trips data.
# 12/3/2020: changing the microtype label names and font sizes for the plots
# 1/12/2021: re-organizing file structure to send to Mona 
##################################################
#mywd <- "/Users/lingjin/Dropbox/Research/FHWA_GeoType/Work/Rscripts"
# mywd = "/Users/xiaodanxu/Library/CloudStorage/GoogleDrive-arielinseu@gmail.com/My Drive/GEMS/typology/TransportGeography-Revision/"
# #mywd <- "~/Box/FHWA"
# setwd(mywd)
# source('initialization.R')
# source('functions.R')

#figuredir <- "./Figures/TransGeo"
#datadir1  <-  "../Data/JGEOTRANS_PaperRevision"
#datadir <- "./Data/CleanData/"
#datadir1 <- "./Data/CleanData/TransGeo"
#rdatadir <-  "./RData/TransGeo"
#tabdir <- './Tables/TransGeo'


# LJ: comment out the following
# # set working directory and sub-directories
# mywd <- "~/Box/FHWA/Task1/TransportGeography-Revision"
# setwd(mywd)
# 
figuredir <- "./Figures"
inputsdir <- "./Data/InputData"
datadir <- "./Data/OutputData"
rdatadir <-  "./Data/RData"
tabdir <- './Tables'
resultsdir <- './Results'
# 
# # Load prep files
# source('initialization.R')
# #source('functions.R')

######
# LOAD DATA
##########
# Compiled inputs for microtypes
FHWA_data <- fread(file.path(inputsdir,"microtypes_inputs.csv")) %>%
  mutate(tract = as.character(tract))

# Commute distance inputs
# Note this is used just for summary stats
#commutes <- read.csv(file.path(inputsdir,"commute_distances_euclidean_tract_updated.csv")) %>% 
 # select(dist_bin, pct_trips_bin) %>%
#  split(. $dist_bin) %>% 
 # purrr::walk(~ stargazer(., type = "text"))

#commutes %>%
 # group_by(dist_bin) %>% 
  #summarise(mean = mean(pct_trips_bin)) %>%
  #stargazer(commutes, type = 'text') 

# load additional inputs
load(file.path(rdatadir,'additional.inputs.RData'))

add.input = add.input %>%
  mutate(tract = as.character(tract))

# Crosswalk to CBSAs and counties
xwalk <- fread(file.path(inputsdir, 'us_xwalk_tract_2017_withID.csv')) %>%
  mutate(tract = as.character(tract))

##################
## # PREPPING DATA FOR CLUSTERING
FHWA_data = FHWA_data %>%
  left_join(add.input) %>%
  select(-w_tract_total, -h_tract_total)

rownames(FHWA_data) <- FHWA_data$tract # set tract as row name for identifier
na.count <- colSums(is.na(FHWA_data)) # check number of missing values for each variables
max(na.count/dim(FHWA_data)[1]) # max missing is about 1.1%, not too bad.

# for now exclude missing values from clustering 
#(we could first derive cluster centroid based on tracts with no missing data on pct trips by distance bins, then assign the tracts based on observed dimensions)
# removed 587 tracts. 0.8% of the tracts.
#NOTE: Sep 10 2021, exporting a data set to impute clusters for the microtypes 

export = as.data.frame(FHWA_data) %>%
  #filter(water == 0) %>%
  dplyr::select(-c(st_code, fips_state, fips_county, cbsa, water)) # remove identifiers 
 
# standardize numeric variables, center and scale them
for (i in 1:length(colnames(export))){
  if (is.numeric(export[1, i])){
    export[, i] <- as.numeric(scale(export[, i]))
  }else{
    export[, i] <- export[, i]
  }
}

write.csv(export, file.path(datadir, "microtypes_inputs_trans_geo_scaled.csv"), row.names = F)

############
# From here down, code is that same as the Task 1 report
FHWA_data = FHWA_data[complete.cases(FHWA_data[,c("pct_trips_bin1",
                                                  "pct_trips_bin2" ,        
                                                 "pct_trips_bin3" ,
                                                 "pct_trips_bin4")])]
names(FHWA_data)
id_vars <- FHWA_data[,c(1:5)] 

# NOTE: dplyr and MASS interfere with the select function
# remove tracts that are 100% water
data_unscaled <- FHWA_data %>%
  filter(water == 0) %>%
  dplyr::select(-c(st_code, fips_state, fips_county, cbsa, water)) # remove identifiers 

rownames(data_unscaled) <- data_unscaled$tract 

data_unscaled = data_unscaled %>%
  select(-tract)

# standardize numeric variables, center and scale them
data_scaled <- as.data.frame(data_unscaled)
for (i in 1:length(colnames(data_scaled))){
  if (is.numeric(data_scaled[1, i])){
    data_scaled[, i] <- as.numeric(scale(data_scaled[, i]))
  }else{
    data_scaled[, i] <- data_scaled[, i]
  }
}
names(data_scaled)

# impute missing data using medians 
data_scaled <- imputeMissings(data_scaled)
raw.df <- data_scaled # export to match Ling's code
save(raw.df, file = file.path(rdatadir, 'raw.scaled.RData'))


#####################
# TABLE 7: Summary stats for first stage inputs
#####################
stargazer(data_unscaled, omit.summary.stat = c("p25", "p75"),  title = "Descriptive Statistics", 
         out = file.path(tabdir,"Sum_stats_stage1_vars.tex"))

# general correlations between variables at national level 
cor_all <- cor(data_scaled, use = "pairwise.complete.obs") #everything
colnames(cor_all)
# label variables
colnames(cor_all) <- c("Broadband", "Non-attainment pollutants", "Land use - Ag", "Pct water",
                       "Development Intensity", "Avg. Circuity", "Dead-end Proportion", "Intersection Density", 
                       "Self-loop Proportion", "Street Density", "Avg Street Length", "Pct manufacturing jobs", 
                       "Pct mining jobs", "Avg IRI", "Pct Full Access Control", "Pct Partial Access Control", 
                       "Pct Truck AADT", "Truck AADT per Lane Mile",  "Lane-miles per Sq. Km" , "Population Density", 
                       "Job Density", "Road grade",  "Lane-meters per Capita","Jobs-Housing Balance", 
                       "Pct Highways" , "Pct local roads", "Pct Midsize Roads",
                       "Pct Trips <1.3 miles",
                       'Pct Trips 1.3-3 miles',
                       'Pct Trips 3-8 miles',
                       'Pct Trips >8 miles',
                       'Trip Source Magnitude')

rownames(cor_all) <-colnames(cor_all)

#######################
# FIGURE 6: Correlation coefficients for microtype inputs
#########################
corrplot(cor_all, tl.col = "black", tl.cex = .7, order = "AOE", # order them in terms of correlations
         method = "circle", type = "upper") # label size = tl.cex

#####################
## DIMENSIONALITY REDUCTION : EFA
######################

##############
# FIGURE 8: SCREE PLOT
##################
png(file = file.path(figuredir, "scree_plot_microtypes.png"))
#Parallel analysis to select the number of factors
# Parallel analysis suggests that the number of factors =  12  and the number of components =  NA 

parallel <- fa.parallel(data_scaled, fm = 'minres', # factor method = min res (does not assume normal distribution)
                        fa = 'fa',  use = 'pairwise', show.legend = TRUE) # missing data
dev.off()
#########################

#Communalities – This is the proportion of each variable's variance that can be explained by the factors 
#(e.g., the underlying latent continua). It is noted as h2 and can be defined as the sum of squared 
# factor loadings for the variables.

# FACTOR ANALYSIS
#In this case, we use select oblique rotation (“oblimin”) as we believe that there is correlation in the factors. 
# (Varimax rotation is used under the assumption that the factors are completely uncorrelated.)
# `Ordinary Least Squared/Minres` factoring (“minres”), is known to provide results similar 
# `Maximum Likelihood` without assuming multivariate normal distribution and
# derives solutions through iterative eigen decomposition like principal axis.

# 10 factor explained 52% variance, the scree plot suggests 10-13 factors are needed.
# we decided 12 factors are more reasonable

nfactor = 12; cutoff = 0.3
for(nfactor in 10:12){
factor = stage1.factor(data_scaled = data_scaled,nfactor = nfactor)
print(factor$loadings,cutoff = cutoff) # show factor loadings and variation explained
print(factor)
}


# save dataset of factors and loadings. nfactor = 12
data_scaled <- cbind(data_scaled, factor$scores)

# export loadings for each tract for visualization
tract_fac <- data_scaled %>%
  select(matches("MR")) %>%
  mutate(tract = row.names(data_unscaled))

write.csv(tract_fac, file = file.path(rdatadir,'tract_loadings.csv'), row.names = F)

###################
# clustering the factors 
factors <- c("MR1", "MR2", "MR3","MR4", "MR5", "MR6", "MR7", "MR8", "MR9", "MR10" ,"MR11", "MR12")

# add individual variables to factors
cluster.df <- data_scaled[factors] 
save(cluster.df, file = file.path(rdatadir,'factor12.scores.RData'))

# uncomment this to load previous factor analysis output
# load(file.path(rdatadir,'factor12.scores_2020ver.RData'))
####################
# plot the 12 factors
factor_names = colnames(factor$loadings)
#[1] "MR2"  "MR1"  "MR4"  "MR5"  "MR10" "MR7"  "MR8"  "MR6"  "MR3"  "MR12" "MR11" "MR9" 
# 060623:  "MR1"  "MR2"  "MR4"  "MR10" "MR5"  "MR8"  "MR12" "MR7"  "MR11" "MR3"  "MR6"  "MR9" 
# TBD
factor_labels = c(
                  'Network Density', # 1
                  'Highway', # 2
                  'Walk/bike Potential', #4
                  'Freight', #10
                  'Jobs Opportunity', # 5
                  'Self-loops', # 8 
                  'Pollution IRI', # 12
                  'Median Trips', #7 
                  'Long Street', # 11
                  'Circuity', # 3
                  'Local Roads', # 6
                  'Job Density' #9
                  ) 

pdf(file = file.path(figuredir,'factor12loadings.cut0.3.pdf'), height = 8, width = 6)
plot.loadings2(factor, factor_labels, names_ordered = transgeo_names_ordered)
dev.off()

#############
# CLUSTERING
################

# Select number of clusters (refer to code: transgeo_stage1_clusterValidation.R)
# results suggest between 6 clusters 
set.seed(1) # so results are replicable each time
tic()
cluster6 <- clara(cluster.df, 6, metric = "euclidean", #euclidean distance metric
                  stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
toc()



# attaching clusters back to raw data
data_scaled <- data_scaled %>% 
  mutate(cluster6 = as.factor(cluster6$clustering),
         # adding identifier back in (do this after reshape)
         tract = row.names(data_unscaled)) %>%
  left_join(xwalk) # link back to county/cbsa etc location info xwalk

write.csv(data_scaled, file = file.path(datadir,'clustering_outputs.csv'), row.names = F)

#####################
# FIGURE 13:  Spider/radar/polar charts of factors by microtype cluster
###########################
spider <- c(factor_names)

# make dataframe of mean of each variable for each cluster (remove all factor variables)
# keep cluster6 variable as factor
radar_data <-  aggregate(data_scaled[spider], list(data_scaled$cluster6), mean)

#radar_data <- aggregate(raw_variables, list(data_scaled$cluster6), mean )
radar_data$Group.1 <- as.numeric(radar_data$Group.1)
colnames(radar_data) <- c("Microtype",factor_labels)

# To use fmsb package, have to add 2 lines to the dataframe: the max and min of each variable to show on the plot
colMax <- function(radar_data) sapply(radar_data, max, na.rm = TRUE) #max
colMin <- function(radar_data) sapply(radar_data, min, na.rm = TRUE) #min

# first row is the column max, second row is the column min
radar_data <- rbind(colMax(radar_data), colMin(radar_data), radar_data)
rownames(radar_data) = c("max", "min", "1", "2", "3", "4", "5", "6")

#colors: darkest to lighest
# 1-#bd0026, 2- #f03b20, 3-#fd8d3c, 4-#feb24c, 5-#fed976, 6- #ffffb2
# make sure centerzero is set correctly, should be the minimum value observed
# KEY: CLUSTER-MICROTYPE: 

png(file = file.path(figuredir,"spider_12fac_clara6_cluster1.png"),height = 700, width = 700)
spider1 <- radarchart(radar_data[1:3, -c(1)], axistype = 0, seg = 4, pty = 32, pdensity = NULL, pangle = 45, 
                     pcol = "#636363", pfcol="#feb24c", plwd=4 , cex.main = 2, #title font size
                      cglty = 3, cglwd = 1, title = "Cluster 1", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, 
                      vlabels = NULL, vlcex = 1.2, caxislabels = NULL, calcex = NULL, paxislabels = NULL, palcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster2.png"),height = 700, width = 700)
spider2 <- radarchart(radar_data[c(1:2,4), -c(1)], axistype = 0, pty = 32, cex.main = 2, pcol = "#636363", pfcol="#f03b20" , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Cluster 2", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster3.png"),height = 700, width = 700)
spider3 <- radarchart(radar_data[c(1:2,5), -c(1)], axistype = 0, pty = 32, pcol = "#636363", cex.main = 2, pfcol="#ffffb2" , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Cluster 3", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster4.png"),height = 700, width = 700)
spider4 <- radarchart(radar_data[c(1:2,6), -c(1)], axistype = 0, pty = 32, pcol = "#636363", cex.main = 2, pfcol="#fd8d3c"  , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Cluster 4", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster5.png"),height = 700, width = 700)
spider5 <- radarchart(radar_data[c(1:2,7), -c(1)], axistype = 0,  pty = 32, pcol = "#636363", cex.main = 2, pfcol="#fed976"  , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Cluster 5", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

png(file = file.path(figuredir,"spider_12fac_clara6_cluster6.png"),height = 700, width = 700)
spider6 <- radarchart(radar_data[c(1:2,8), -c(1)], axistype = 0,  pty = 32, pcol = "#636363", cex.main = 2, pfcol="#bd0026"  , plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "Cluster 6", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
dev.off()

# plotting them all together:
colors_reds <- c("#feb24c",  "#f03b20", "#ffffb2", "#fd8d3c","#fed976", "#bd0026")
png(file = "./Figures/spider_10fac_clara6_all.png")
spider_all <- radarchart(radar_data[, -c(1)], axistype = 0,  pty = 32, plty = 1, pcol = colors_reds, plwd=4 , palcex = 1.2,
                      cglty = 3, cglwd = 1, title = "All microtypes", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.4)
legend(x=1, y=1, legend = rownames(radar_data[-c(1,2),]), bty = "n", pch=20 , col=colors_reds , text.col = "black", cex=1.2, pt.cex=3)
dev.off()

#export microtypes for use in second-stage clusters

# add cluster to ID datasets
data_scaled$tract <- row.names(data_unscaled)
keep = data_scaled %>% 
  select(tract, cluster6)

plot_data = id_vars %>%
  left_join(keep) %>%
  mutate(cluster6 = as.numeric(cluster6))

# add type 7 for tracts that have only water
plot_data$cluster6[is.na(plot_data$cluster6)] <- 7

#export data to use for geotypes development
write.csv(plot_data, file = file.path(datadir,"microtypes_output_clara6.csv"), row.names = F)

#################
# TABLE 9: mean inputs by microtpye 
###################
results <- plot_data %>% 
  select(tract, cluster6) %>% 
  filter(cluster6 <7)

sum_stats = FHWA_data %>%
  merge(results, by = "tract")

sum_stats6 <- sum_stats %>% 
  group_by(cluster6) %>% 
  select(-tract,-st_code, -fips_county, - cbsa, -water) %>% 
  summarise_all("mean", na.rm = T)

write.csv(sum_stats6, file = file.path(tabdir, "sumstats_by_microtype_clara6.csv"), row.names = F)




