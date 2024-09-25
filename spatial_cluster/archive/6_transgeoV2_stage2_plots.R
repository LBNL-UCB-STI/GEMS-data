# Create box plots by geotype

#####
# LOAD DATA
#######
results <- read.csv(file.path(resultsdir, "stage2_clustering_outputs_pam6.csv"))

# Crosswalk between tracts and counties/CBSA
xwalk <- read.csv(file.path(inputsdir, "us_xwalk_tract_2017_withID.csv"))
###############


means <- results %>% 
  group_by(cluster6) %>%
  summarise(mean_size = mean(fhwa_type_num))

stargazer(results, type = 'text')

#keep only one value per spatial ID
results <- results %>% 
  select(spatial_id, fhwa_type_num, pct_micro.6, hhi_normalized, emp_centers, cluster6) %>%
  distinct() %>%
  filter(!is.na(cluster6) & !is.na(hhi_normalized)) %>%
  mutate(emp_scaled = scale(emp_centers),
         emp_norm = (emp_centers-min(emp_centers))/(max(emp_centers)-min(emp_centers)),
         geotype = case_when(
            cluster6 == 6 ~ 'A' ,
            cluster6 == 5 ~ 'B',
            cluster6 == 2 ~ 'C',
            cluster6 == 3 ~ 'D',
            cluster6 == 1 ~ 'E', 
            cluster6 == 4 ~ 'F'       ))

#basic sum stats
stargazer(results, type = 'text')


#export labeled file of all results 

#cbsa results
cbsa <- results %>%
  merge(xwalk, by = "spatial_id") %>%
  filter(cbsa <99999) %>%
  select(spatial_id, cbsaname, cbsa, pct_micro.6, hhi_normalized, emp_centers, geotype) %>%
  distinct()

table(cbsa$geotype)
write.csv(cbsa, file.path(resultsdir, "cbsa_results_labeled.csv"), row.names = F)

#rural 
rural <- results %>%
  merge(xwalk, by = "spatial_id") %>%
  filter(cbsa ==99999) %>%
  select(spatial_id, pct_micro.6, hhi_normalized, emp_centers, geotype) %>%
  distinct()

write.csv(rural, file.path(resultsdir, "rural_results_labeled.csv"), row.names = F)


box <- c("pct_micro.6", "hhi_normalized", "emp_norm", "geotype")

box.df <- reshape2::melt(results[box], id.var = "geotype") # no demographic

# GEOTYPE colors: darkest to lightest
geo_colors <- c("#014636","#016c59", "#02818a", "#3690c0", "#67a9cf", "#a6bddb")
# A - #014636, 
# B - #016c59,
# C - #02818a, 
# D - #3690c0, 
#E - #67a9cf, 
#F - #a6bddb, 

# faceted box plots by variables
labels <- c("Proportion of Microtype 1", "Commute Dispersion", "Polycentricity")


png(file.path(figuredir,"geo_boxplot_trans_geo.png"), height = 400, width = 700)
ggplot(data = box.df, aes(x=variable, y=value)) +
  geom_violin(aes(fill=geotype)) +
  scale_fill_manual(values = geo_colors, guide = guide_legend()) +
  guides(fill = guide_legend(legend.position = "bottom",
                             legend.box = "horizontal", 
                             nrow = 1)) +
  theme(legend.position = "bottom", text = element_text(size = 20)) +
  xlab("Input") +
  scale_x_discrete(labels = labels) +
  labs(title="Geo-econ variables, 7 clusters with CLARA") 
dev.off()

# combined boxplot
ggplot(data = box.df, aes(x=cluster6, y=value)) +
  geom_boxplot(aes(fill=variable), outlier.colour=NA) +
#  coord_cartesian(ylim=c(0,50)) +
  scale_fill_brewer(palette = "YlGnBu") +
  theme(legend.position = "bottom", axis.text.x =element_text(size = 10) ) +
  labs(title="Geo-econ variables, 7 clusters with CLARA")
