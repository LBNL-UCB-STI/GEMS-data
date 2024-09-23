# CODE FOR SECOND-STAGE CLUSTERING TO GENERATE MICROTYPES
# IMPORT CLEANED DATASET OF INPUTS
# GENERATE PRINCIPAL COMPONENTS WITH PCA
# CLUSTER USING K-MEANS
# SAVE SHAPEFILES OF GEOOTYPES

# NATALIE POPOVICH
# LAWRENCE BERKELEY NATIONAL LAB
# LAST UPDATED: NOV 6 2020
# LJ 11/30/2020: updated with PAM according to cluster validation
##################################################
inputs <- read.csv(file.path(inputsdir,"geotypes_inputs_transp_geo.csv")) 

# link back to county/cbsa etc location info xwalk
xwalk <- read.csv(file = file.path(inputsdir, "us_xwalk_tract_2017_withID.csv")) %>%
  select(cbsa, cbsaname, spatial_id, st_code, cty, ctyname) %>%
  distinct() 

## # PREPPING DATA FOR CLUSTERING
rownames(inputs) <- inputs$spatial_id # set tract as row name for identifier
inputs$emp_centers <- as.numeric(inputs$emp_centers)
# check number of missing values for each variables, 
# HHI_normalized had 44 missing, which is about 44/2250 = 2%
colSums(is.na(inputs)) 

# take a look at the 44 rows, 
# these are either high% in microtype 2 (the high circuity ones) or all in water
# 0% in densest microtype6.
# Actually, hhi_normalized for all water geotypes are not defined, so clustering exclude the all water geotypes.
# which will exclude 20 counties with only pct_micro.na = 1, and other variables = 0 or not defined.
inputs %>% filter(is.na(hhi_normalized))


# look at rows with high mirotype2, we can see the median hhi_normalized is lower in these  than the pop
#tmp = inputs %>% filter(pct_micro.2>0.5)

# remove all water rows
inputs = inputs %>% filter(pct_micro.na < 1)
colSums(is.na(inputs)) 


# standardize numeric variables, center and scale them
data_unscaled <- inputs[,-c(1)] %>%
  select(-fhwa_type, - fhwa_type_num)
# first check if any column is not numeric
sum(apply(data_unscaled,c(2), is.numeric)) # 9, all numeric

data_scaled <- as.data.frame(scale(data_unscaled)) # scale to mean = 0, sd = 1, returns a matrix, 

# and turn back to data.frame for later imputation function
names(data_scaled)

## get rid of colinear vars before running
# data_scaled <- data_scaled %>% dplyr::select(-matches('7'))
#data_scaled <- data_scaled %>% dplyr::select(-pct_within, -pctdist5_10_total)

##################### clustering (without PCA as there are now very few variables) ######################
# because the hhi in high microtype2 geotypes appeared to be lower than the population average or median, 
# here we use knn to impute so that the missing rows have most similar characteristics across all other covariates.
# k ~ sqrt(2200) ~ 47

rownames(data_scaled) = inputs$spatial_id

train.dat = data_scaled %>% 
  filter(!is.na(hhi_normalized)) %>%
  select(-hhi_normalized)
train.out = data_scaled %>% 
  filter(!is.na(hhi_normalized)) %>%
  select(hhi_normalized)

reg_pred = data_scaled %>% 
  filter(is.na(hhi_normalized)) %>%
  select(-hhi_normalized)

knn_res <- knn.reg(train.dat, reg_pred, train.out, k = 47)

data_scaled[is.na(data_scaled$hhi_normalized),'hhi_normalized'] = knn_res$pred

save(data_scaled, file = file.path(rdatadir,"raw_scaled_stage2.RData"))


# general correlations between all geotype inputs
cor_all <- cor(data_scaled, use = "pairwise.complete.obs") #everything
colnames(cor_all)
# label variables
# check the order of these variables! will need to re-label the clusters with new microtypes
colnames(cor_all) <- c("Pct. cluster 1", "Pct. cluster 2", 
                       "Pct. cluster 3", "Pct. cluster 4", 
                       "Pct. cluster 5", "Pct. cluster 6",
                       "Pct. water",
                       "HHI of commutes",  "No. of employment centers")

#  c("Pct. Microtype 1", "Pct. Microtype 2", "Pct. Microtype 3", "Pct. Microtype 6", "Pct. Microtype 4", "Pct. Microtype 5",
#                       "HHI of commutes",  "No. of employment centers")
rownames(cor_all) <-colnames(cor_all)

##############3
# FIGURE 7: Correlation coefficients between geotype inputs 
###############
corrplot(cor_all, tl.col = "black", tl.cex = .7, order = "AOE", # order them in terms of correlations
         method = "circle", type = "upper") # label size = tl.cex

###############
# TABLE 8: Summary statistics of second-stage inputs
#################
stargazer(inputs[, -c(1)], omit.summary.stat = c("p25", "p75"),  title = "Descriptive Statistics", 
          out = file.path(tabdir,"Sum_stats_regional.tex"))

##############
# No DIMENSIONALITY REDUCTION
##################
#################
# CONDUCT CLUSTERING 
#####################
# First determine the number of clusters


testK.kmeans = testK(cluster.df = data_scaled%>%select(pct_micro.6,hhi_normalized,emp_centers),
                     k.range = 2:10, method = 'kmeans')

testK.pam = testK(cluster.df = data_scaled%>%select(pct_micro.6,hhi_normalized,emp_centers), 
                  k.range = 2:10, method = 'pam')

ggplot(testK.kmeans%>%filter(k<=8), aes(x = k, y= davies_bouldin)) +
  theme_bw()+
  geom_vline(xintercept = 5, size = rel(1.1), color = 'grey') +
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
  scale_x_continuous(breaks = 2:8, labels = 2:8)+
  ggtitle('Method = K-Means')

ggsave(file = file.path(figuredir,'stage2_cluster_kmeans_validity_pick5.pdf'),
       height = 3, width = 5) 


ggplot(testK.pam%>%filter(k<=8), aes(x = k, y= davies_bouldin)) +
  theme_bw()+
  geom_vline(xintercept = 6, size = rel(1.1), color = 'grey') +
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
  scale_x_continuous(breaks = 2:8, labels = 2:8)+
  ggtitle('Method = PAM')


ggsave(file = file.path(figuredir,'stage2_cluster_pam_validity_pick6.pdf'),
       height = 3, width = 5) 


### Now pick k = 6 and run 
set.seed(10)
k=6
#cluster <- pam(data_scaled, 6) 

cluster <- pam(data_scaled%>%select(pct_micro.6,hhi_normalized, emp_centers), k) 

inputs$cluster6 <- as.factor(cluster$clustering)

table(inputs$cluster6)

table(inputs$cluster6, inputs$fhwa_type_num)



# merge cluster results with the spatial IDs for all counties
# NOTE: update this with the new number of geotypes

clusters <- xwalk %>%
  left_join(inputs, by = "spatial_id")

write.csv(clusters, file = file.path(resultsdir,'stage2_clustering_outputs.csv'), row.names = F)

#export csv of results separately by county and cbsa
cbsa <- xwalk %>% 
  select(cbsaname, spatial_id, cbsa) %>% 
  filter(cbsa < 99999) %>%
  distinct()

cbsa_results <- cbsa %>%
  left_join(inputs, by = "spatial_id")

write.csv(cbsa_results, file.path(resultsdir,'geotypes_cbsa_pam6.csv'), row.names = F)

rural <- xwalk %>% 
  select(spatial_id, cty, ctyname, st_code, cbsa) %>% 
  filter(cbsa == 99999)

rural_results <- rural %>%
  left_join(clusters, by = "spatial_id")

write.csv(rural_results, file.path(resultsdir,'geotypes_rural_pam6.csv'), row.names = F)

##################
# Table 10: median input values by geotype (Not in paper)
########################
inputs <- inputs %>% 
  select(-spatial_id)

sum_stats <- inputs %>% 
  select('pct_micro.6','hhi_normalized','emp_centers',cluster6)%>%
  group_by(cluster6) %>% 
  summarise_all(c("IQR","median"), na.rm = T)

write.csv(sum_stats, file = file.path(tabdir,"sumstats_by_geotype_pam6.csv"), row.names = F)
