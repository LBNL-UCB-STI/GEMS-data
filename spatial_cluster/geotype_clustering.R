#inputsdir <- "/Users/jessicalazarus/Documents/GEMS/Input Data"

mywd <- "C:/Users/xiaodanxu/Documents/GitHub/GEMS-data/spatial_cluster"
setwd(mywd)
source('initialization.R')
source('functions.R')

# library(cluster)
# library(tidyverse)
# library(data.table)
# library(dplyr)
# library(tidyr)

demandinputsdir <- "C:/FHWA_R2/Demand/CleanData"
spatialinputsdir <- "C:/FHWA_R2/spatial_boundary/CleanData"
figuredir <- "C:/FHWA_R2/Demand/Figures"
outputsdir <- "C:/FHWA_R2/Network/CleanData"

inputs <- read.csv(file.path(demandinputsdir,"geotype_inputs.csv"))

# link back to county/cbsa etc location info xwalk
xwalk <- read.csv(file = file.path(spatialinputsdir, "cleaned_lodes8_crosswalk_with_ID.csv")) %>%
  select(cbsa, cbsaname, spatial_id, st, cty, ctyname) %>%
  distinct() 

## # PREPPING DATA FOR CLUSTERING
rownames(inputs) <- inputs$spatial_id # set tract as row name for identifier

# combine correlated variables
inputs$office_service_jobs <- inputs$office_jobs + inputs$service_jobs
inputs$recreation_retail_jobs <- inputs$retail_jobs + inputs$recreation_jobs

## get rid of colinear vars before running
inputs <- inputs %>% dplyr::select(-office_jobs, -retail_jobs, -recreation_jobs, -service_jobs, -laneMiles, -ALAND, -AWATER, -total_jobs )


# check number of missing values for each variables, 
colSums(is.na(inputs)) 
inputs = inputs %>% filter(!is.na(recreation_retail_jobs))
colSums(is.na(inputs))

# standardize numeric variables, center and scale them
# first check if any column is not numeric
sum(apply(inputs,c(2), is.numeric)) # all numeric
inputs_cbsa = inputs %>% filter(is_cbsa==1) %>% dplyr::select(-is_cbsa)
inputs_ncbsa = inputs %>% filter(is_cbsa==0) %>% dplyr::select(-is_cbsa)
data_scaled_cbsa <- as.data.frame(scale(inputs_cbsa)) %>% dplyr::select(-spatial_id)# scale to mean = 0, sd = 1, returns a matrix,
data_scaled_ncbsa <- as.data.frame(scale(inputs_ncbsa)) %>% dplyr::select(-spatial_id)# scale to mean = 0, sd = 1, returns a matrix,

# and turn back to data.frame for later imputation function
names(data_scaled_cbsa)
rownames(data_scaled_cbsa) = inputs_cbsa$spatial_id
rownames(data_scaled_ncbsa) = inputs_ncbsa$spatial_id
sum(apply(data_scaled_cbsa,c(2), is.numeric))
##################### clustering (without PCA as there are now very few variables) ######################
# general correlations between all geotype inputs
cor_all <- cor(data_scaled_cbsa, use = "pairwise.complete.obs") #everything
colnames(cor_all)
corrplot(cor_all, tl.col = "black", tl.cex = .7, order = "AOE",
         method = "circle", type = "upper")
cor_all <- cor(data_scaled_ncbsa, use = "pairwise.complete.obs") #everything
corrplot(cor_all, tl.col = "black", tl.cex = .7, order = "AOE", # order them in terms of correlations
         method = "circle", type = "upper") # label size = tl.cex

####################
# CONDUCT CLUSTERING 
#####################
# First determine the number of clusters
max_c = 10
testK.kmeans = testK(cluster.df = data_scaled_ncbsa,
                     k.range = 2:max_c, method = 'kmeans')

testK.pam = testK(cluster.df = data_scaled_ncbsa, 
                  k.range = 2:max_c, method = 'pam')

ggplot(testK.kmeans%>%filter(k<=max_c), aes(x = silhouette, y= davies_bouldin)) +
  theme_bw()+
  geom_point() +
  scale_y_continuous(
    # Features of the first axis
    name = "Inverse DBI")+
  theme(
    axis.title.y = element_text( size=13)) +
  xlab('Silhouette') +
  ggtitle('Method = K-Means')
ggsave(file = file.path(figuredir,'stage2_cluster_ncbsa_kmeans_pareto.png'),
       height = 3, width = 5) 

ggplot(testK.kmeans%>%filter(k<=max_c), aes(x = k, y= davies_bouldin)) +
  theme_bw()+
  #geom_vline(xintercept = 5, size = rel(1.1), color = 'grey') +
  geom_point(color = '#f59ff2') +
  geom_line(color = '#f59ff2') +
  geom_point(aes(x = k, y = silhouette*4), color = "#7be2ed")+
  geom_line(aes(x = k, y = silhouette*4),color = "#7be2ed") +
  scale_y_continuous(
    # Features of the first axis
    name = "Inverse DBI",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./4, name="Silhouette")) + 
  theme(
    axis.title.y = element_text(color ='#f59ff2', size=13),
    axis.title.y.right = element_text(color = "#7be2ed", size=13)) +
  xlab('Number of Clusters') +
  scale_x_continuous(breaks = 2:max_c, labels = 2:max_c)+
  ggtitle('Method = K-Means')

ggsave(file = file.path(figuredir,'stage2_cluster_ncbsa_kmeans_validity.png'),
       height = 3, width = 5) 

ggplot(testK.pam%>%filter(k<=max_c), aes(x = silhouette, y= davies_bouldin)) +
  theme_bw()+
  geom_point() +
  scale_y_continuous(
    # Features of the first axis
    name = "Inverse DBI")+
  theme(
    axis.title.y = element_text( size=13)) +
  xlab('Silhouette') +
  ggtitle('Method = K-Means')
ggsave(file = file.path(figuredir,'stage2_cluster_ncbsa_pam_pareto.png'),
       height = 3, width = 5) 

ggplot(testK.pam%>%filter(k<=max_c), aes(x = k, y= davies_bouldin)) +
  theme_bw()+
  #geom_vline(xintercept = 6, size = rel(1.1), color = 'grey') +
  geom_point(color = '#f59ff2') +
  geom_line(color = '#f59ff2') +
  geom_point(aes(x = k, y = silhouette*4), color = "#7be2ed")+
  geom_line(aes(x = k, y = silhouette*4),color = "#7be2ed") +
  scale_y_continuous(
    # Features of the first axis
    name = "Inverse DBI",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./4, name="Silhouette")) + 
  theme(
    axis.title.y = element_text(color ='#f59ff2', size=13),
    axis.title.y.right = element_text(color = "#7be2ed", size=13)) +
  xlab('Number of Clusters') +
  scale_x_continuous(breaks = 2:max_c, labels = 2:max_c)+
  ggtitle('Method = PAM')

ggsave(file = file.path(figuredir,'stage2_cluster_ncbsa_pam_validity.png'),
       height = 3, width = 5) 

### Now pick k and run 
set.seed(10)
k=2

cluster_cbsa <- pam(data_scaled_cbsa, k) 
inputs_cbsa <- inputs_cbsa %>% mutate(cluster2=as.factor(cluster_cbsa$cluster))
cluster_ncbsa <- pam(data_scaled_ncbsa, k) 
inputs_ncbsa <- inputs_ncbsa %>% mutate(cluster2=as.factor(cluster_ncbsa$cluster))
data_scaled_cbsa <- data_scaled_cbsa %>% mutate(cluster2=as.factor(cluster_cbsa$cluster)) 
data_scaled_ncbsa <- data_scaled_ncbsa %>% mutate(cluster2=as.factor(cluster_ncbsa$cluster))

table(inputs_cbsa$cluster2)
table(inputs_ncbsa$cluster2)

k=3

cluster_cbsa <- pam(data_scaled_cbsa, k) 
inputs_cbsa <- inputs_cbsa %>% mutate(cluster3=as.factor(cluster_cbsa$cluster))
cluster_ncbsa <- pam(data_scaled_ncbsa, k) 
inputs_ncbsa <- inputs_ncbsa %>% mutate(cluster3=as.factor(cluster_ncbsa$cluster))
data_scaled_cbsa <- data_scaled_cbsa %>% mutate(cluster3=as.factor(cluster_cbsa$cluster)) 
data_scaled_ncbsa <- data_scaled_ncbsa %>% mutate(cluster3=as.factor(cluster_ncbsa$cluster))

table(inputs_cbsa$cluster3)
table(inputs_ncbsa$cluster3)

k=4

cluster_cbsa <- pam(data_scaled_cbsa, k) 
inputs_cbsa <- inputs_cbsa %>% mutate(cluster4=as.factor(cluster_cbsa$cluster))
cluster_ncbsa <- pam(data_scaled_ncbsa, k) 
inputs_ncbsa <- inputs_ncbsa %>% mutate(cluster4=as.factor(cluster_ncbsa$cluster))
data_scaled_cbsa <- data_scaled_cbsa %>% mutate(cluster4=as.factor(cluster_cbsa$cluster)) 
data_scaled_ncbsa <- data_scaled_ncbsa %>% mutate(cluster4=as.factor(cluster_ncbsa$cluster))

table(inputs_cbsa$cluster4)
table(inputs_ncbsa$cluster4)

write.csv(inputs_cbsa, file.path(outputsdir,"Urban_Clusters_k2_4_pam.csv"))
write.csv(inputs_ncbsa, file.path(outputsdir,"Rural_Clusters_k2_4_pam.csv"))
write.csv(data_scaled_cbsa, file.path(outputsdir,"Urban_Clusters_k2_4_pam_scaled.csv"))
write.csv(data_scaled_ncbsa, file.path(outputsdir,"Rural_Clusters_k2_4_pam_scaled.csv"))


cluster_cbsa <- pam(data_scaled_cbsa, 7) 
inputs_cbsa <- inputs_cbsa %>% mutate(cluster7=as.factor(cluster_cbsa$cluster))
cluster_ncbsa <- pam(data_scaled_ncbsa, 8) 
inputs_ncbsa <- inputs_ncbsa %>% mutate(cluster8=as.factor(cluster_ncbsa$cluster))
data_scaled_cbsa <- data_scaled_cbsa %>% mutate(cluster7=as.factor(cluster_cbsa$cluster)) 
data_scaled_ncbsa <- data_scaled_ncbsa %>% mutate(cluster8=as.factor(cluster_ncbsa$cluster))


table(inputs_cbsa$cluster7)
table(inputs_ncbsa$cluster8)

write.csv(inputs_cbsa, file.path(outputsdir,"Urban_Clusters_k7_pam.csv"))
write.csv(inputs_ncbsa, file.path(outputsdir,"Rural_Clusters_k8_pam.csv"))
write.csv(data_scaled_cbsa, file.path(outputsdir,"Urban_Clusters_k7_pam_scaled.csv"))
write.csv(data_scaled_ncbsa, file.path(outputsdir,"Rural_Clusters_k8_pam_scaled.csv"))
k <- 2
cluster_cbsa <- kmeans(data_scaled_cbsa, k) 
inputs_cbsa <- inputs_cbsa %>% mutate(cluster2=as.factor(cluster_cbsa$cluster))
cluster_ncbsa <- kmeans(data_scaled_ncbsa, k) 
inputs_ncbsa <- inputs_ncbsa %>% mutate(cluster2=as.factor(cluster_ncbsa$cluster))

table(inputs_cbsa$cluster2)
table(inputs_ncbsa$cluster2)

# # Access centroids (optional)
# centroids_cbsa <- cluster_cbsa$centers
# centroids_ncbsa <- cluster_ncbsa$centers
# 
# # Create radar charts for each cluster
# library(fmsb)  # for radar charts
# par(mfrow=c(1,k))  # arrange plots in a single row
# 
# spider <- c(list_of_attr, 'demand_microtype_comb')
# list_attr <- c("pct_water","agriculture_frac","industry_jobs","education_jobs","healthcare_jobs","government_jobs","job_density","HHI","lm_density","Rural_agriculture" ,"Rural_towncenter","Suburban","Urban_center","Urban_industrial","pop_density","office_service_jobs","recreation_retail_jobs")
# # make dataframe of mean of each variable for each cluster (remove all factor variables)
# # keep cluster6 variable as factor
# radar_data <- inputs_cbsa
# radar_data <-  aggregate(radar_data[list_attr], list(inputs_cbsa$cluster2), median)
# 
# #radar_data <- aggregate(raw_variables, list(data_scaled$cluster6), mean )
# #radar_data$Group.1 <- as.numeric(radar_data$Group.1)
# #colnames(radar_data) <- c("Microtype",factor_labels)
# 
# # To use fmsb package, have to add 2 lines to the dataframe: the max and min of each variable to show on the plot
# colMax <- function(radar_data) sapply(radar_data, max, na.rm = TRUE) #max
# colMin <- function(radar_data) sapply(radar_data, min, na.rm = TRUE) #min
# 
# # first row is the column max, second row is the column min
# #radar_data <- rbind(colMax(radar_data), colMin(radar_data), radar_data)
# #rownames(radar_data) = c("max", "min", "Rural_agriculture", "Rural_towncenter", "Suburban", "Urban_center", "Urban_industrial")
# 
# radar_data <- radar_data  %>% select(-Group.1 )
# radar_data  <- mutate_all(radar_data , function(x) as.numeric(as.character(x)))
# 
# png(file = file.path(figuredir,"spider_urban_2clusters_cluster1.png"),height = 700, width = 700)
# #spider1 <- radarchart(radar_data, axistype = 0, #seg = 4, pty = 32, pdensity = NULL, pangle = 45, 
# #                      pcol = "#636363", pfcol="#0000B3FF", plwd=4 , cex.main = 2, #title font size
# #                      cglty = 3, cglwd = 1, title = "Urban 1", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, 
# #                      vlabels = NULL, vlcex = 1.2, caxislabels = NULL, calcex = NULL, paxislabels = NULL, palcex = 1.2)
# spider1 <- radarchart(radar_data, axistype = 0, #seg = 4, pty = 32, pdensity = NULL, pangle = 45, 
#                       title = "Urban 1"
#                       )
# 
# dev.off()
# 
# png(file = file.path(figuredir,"spider_12fac_clara6_cluster2.png"),height = 700, width = 700)
# spider2 <- radarchart(radar_data[c(1:2,4), ], axistype = 0, pty = 32, cex.main = 2, pcol = "#636363", pfcol="#4500FFFF" , plwd=4 , palcex = 1.2,
#                       cglty = 3, cglwd = 1, title = "Rural towncenter", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
# dev.off()
# 
# png(file = file.path(figuredir,"spider_12fac_clara6_cluster3.png"),height = 700, width = 700)
# spider3 <- radarchart(radar_data[c(1:2,5), ], axistype = 0, pty = 32, pcol = "#636363", cex.main = 2, pfcol="#C527D8FF" , plwd=4 , palcex = 1.2,
#                       cglty = 3, cglwd = 1, title = "Suburban", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
# dev.off()
# 
# png(file = file.path(figuredir,"spider_12fac_clara6_cluster4.png"),height = 700, width = 700)
# spider4 <- radarchart(radar_data[c(1:2,6), ], axistype = 0, pty = 32, pcol = "#636363", cex.main = 2, pfcol= "#FFA35CFF", plwd=4 , palcex = 1.2,
#                       cglty = 3, cglwd = 1, title = "Urban center", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
# dev.off()
# 
# png(file = file.path(figuredir,"spider_12fac_clara6_cluster5.png"),height = 700, width = 700)
# spider5 <- radarchart(radar_data[c(1:2,7), ], axistype = 0,  pty = 32, pcol = "#636363", cex.main = 2, pfcol=  "#FFF50AFF", plwd=4 , palcex = 1.2,
#                       cglty = 3, cglwd = 1, title = "Urban industrial", maxmin = TRUE, na.itp = TRUE, centerzero = FALSE, vlcex = 1.2)
# dev.off()

