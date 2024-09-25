# LJ 8/20/2024
# summarise the MFD parameters by typology types
# and some visualization

# link to repo
mywd = "C:/Users/xiaodanxu/Documents/GitHub/GEMS-data/mfd"
setwd(mywd)
#source('initialization.R')
#source('functions.R')
library(pacman)
p_load(R.utils) # for reading gz files
p_load(readr) # for read files from a zip foler

p_load(data.table)
p_load(dplyr)
p_load(tidyr)
p_load(ggplot2)
p_load(stringr)
p_load(vroom)
p_load(tictoc)
p_load(parallel)
p_load(doParallel)
p_load(foreach)
p_load(stargazer)
p_load(Hmisc) # weighted quantile



#source('initialization.R')
source('functions.R')
auxdatdir <- 'C:/data/CATTLab_delivery/AuxiliaryData' # the blob-microtype ids
rdatadir <-  "C:/data/RData/MFD_v2"
figuredir <- "C:/data/Figures/MFD_v2"
#outdir  <-  "../Data/"

part.dir <- 'C:/data/CATTLab_delivery/National/proc_data/partition_input'
cat.proc = 'C:/data/CATTLab_delivery/National/proc_data/all_tracts_frc3_5'  # tract level data aggregating frc3-5 consistently
cat.proc2 <- 'C:/data/CATTLab_delivery/National/proc_data/m3_tracts_frc1_5' # only tract level data for microtype 3 aggregating frc1-5

linker = fread(file.path(auxdatdir,'ccst_geoid_key_tranps_geo_with_imputation.csv'))%>% # old mcirotype xwalk
  select(GEOID, microtype, st_code, cbsa, cbsaname,geotype)

###########################
# 1. load network typlogy
ntype <- fread(file = file.path(auxdatdir, 'microtype_geotype_output_2010.csv')) %>%
  mutate(tract_geoid = as.numeric(GEOID))%>%
  select(-cbsa, -cbsaname)

# rename per Anna's requestt
ntype <- ntype %>%
  mutate(geotype = recode(geotype, 
                          'CBSA_2' = 'A', 
                          'CBSA_1' = 'B',
                          'NONCBSA_1' = 'C',
                          'NONCBSA_2' = 'D'),
         network_microtype = recode(network_microtype,
                                    'Urban_1' = 'Urban_5',
                                    'Urban_2' = 'Urban_4',
                                    'Urban_4' = 'Urban_2',
                                    'Urban_5' = 'Urban_1',
                                    'Rural_1' = 'Rural_3',
                                    'Rural_3' = 'Rural_1'))




###########################
# 2. load partition MFD estimates and merge to tract_geoid
load(file = file.path(cat.proc, "partition.mfd.fit.RData")) # mfd.res_steplinear, mfd.res_parabolic, mfd.res_spddens, 
load(file = file.path(part.dir, 'reagg_partition_xwalk.RData')) # load pt
mfd.res_parabolic <- mfd.res_parabolic %>%
  left_join(pt) %>%
  mutate(tract_geoid = as.numeric(GEOID))

# mfd.res_parabolic <- mfd.res_parabolic %>%
#   select()

# 3. load the tract level MFD estimates and exclude the tracts in partitions
load(file = file.path(cat.proc, 'tract.mfd.fit.res.RData')) #fit.rural.all, fit.parabolic.all, fit.spddens.all, 

tract.mfd <- fit.parabolic.all %>%
  mutate(tract_geoid = as.numeric(tract_geoid))%>%
  filter(!(tract_geoid %in% mfd.res_parabolic$tract_geoid))


# 3.2 load the tract level MFD estiamtes for only old microtype 3
load(file = file.path(cat.proc2, 'tract.m3.mfd.fit.res.RData')) #fitm33.rural.all, fitm3.parabolic.all, fitm3.spddens.all, 

tract.m3.mfd <- fitm3.parabolic.all %>%
  mutate(tract_geoid = as.numeric(tract_geoid))%>%
  inner_join(ntype) %>%
  filter(network_microtype == 'Urban_2') # change to Urban_2, which is the new name for Urban_4
# 4. combine the data and summarise for non-freeway tracts, 
# i.e. excluding tracts asssigned to 'urban_4'

mfd.res_parabolic <- mfd.res_parabolic %>%
  ungroup()%>%
  select(names(tract.mfd))

mfd.fit.parabolic <- bind_rows(mfd.res_parabolic, tract.mfd) %>%
  left_join(ntype) %>%
  filter(network_microtype != 'Urban_2')

# bind rows with urban_4 tracts that has data for frc1-5 aggregation
mfd.fit.parabolic = bind_rows(mfd.fit.parabolic,tract.m3.mfd)

median.fit.parabolic <- mfd.fit.parabolic %>%
  group_by(geotype, network_microtype) %>%
  summarise(capacity_plph = median(capcity.fit.parabolic,na.rm = T),
            max_spd_mph = median(max.spd.parabolic, na.rm = T),
            a = - max_spd_mph^2/capacity_plph/4,
            b =  4*capacity_plph/max_spd_mph,
            critical_density_plpm = b/2) 

write.csv(median.fit.parabolic,file = file.path(rdatadir,'median_mfd_parameters_3x5.csv'),row.names = F)

# 5. some visualization to see if parameters vary by geotypes or not
ggplot(mfd.fit.parabolic %>% filter(geotype %in% c('A','B')), aes(x = network_microtype, y = critical.density.fit.parabolic, color = network_microtype))+
  geom_violin(width = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(c(0,100))+
  ylab('critical_density (#veh/lane/mile)')+
  facet_grid(~geotype)+
  ggtitle('CBSA critical_density from parabolic fit')
ggsave(file.path(figuredir, 'urban_ntype.cridensity.parabolic.by.cluster.pdf'), height =3, width = 8.5)

ggplot(mfd.fit.parabolic %>% filter(geotype %in% c('A','B')), aes(x = network_microtype, y = capcity.fit.parabolic, color = network_microtype))+
  geom_violin(width = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(c(0,2500))+
  ylab('max capacity (#veh/lane/hour)')+
  facet_grid(~geotype)+
  ggtitle('CBSA max capacity from parabolic fit')
ggsave(file.path(figuredir, 'urban_ntype.maxcapacity.parabolic.by.cluster.pdf'), height =3, width = 8.5)

ggplot(mfd.fit.parabolic %>% filter(geotype %in% c('A','B')), aes(x = network_microtype, y = max.spd.parabolic, color = network_microtype))+
  geom_violin(width = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(c(0,70))+
  ylab('max speed (mile/hour)')+
  facet_grid(~geotype)+
  ggtitle('CBSA max speed from parabolic fit')
ggsave(file.path(figuredir, 'urban_ntype.maxspd.parabolic.by.cluster.pdf'), height =3, width = 8.5)


## plot rural
ggplot(mfd.fit.parabolic %>% filter(geotype %in% c('C','D')), aes(x = network_microtype, y = critical.density.fit.parabolic, color = network_microtype))+
  geom_violin(width = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(c(0,100))+
  ylab('critical_density (#veh/lane/mile)')+
  facet_grid(~geotype)+
  ggtitle('NONCBSA critical_density from parabolic fit')
ggsave(file.path(figuredir, 'rural_ntype.cridensity.parabolic.by.cluster.pdf'), height =3, width = 7.5)

ggplot(mfd.fit.parabolic %>% filter(geotype %in% c('C','D')), aes(x = network_microtype, y = capcity.fit.parabolic, color = network_microtype))+
  geom_violin(width = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(c(0,2500))+
  ylab('max capacity (#veh/lane/hour)')+
  facet_grid(~geotype)+
  ggtitle('NONCBSA max capacity from parabolic fit')
ggsave(file.path(figuredir, 'rural_ntype.maxcapacity.parabolic.by.cluster.pdf'), height =3, width = 7.5)

ggplot(mfd.fit.parabolic %>% filter(geotype %in% c('C','D')), aes(x = network_microtype, y = max.spd.parabolic, color = network_microtype))+
  geom_violin(width = 1) +
  geom_boxplot(width = 0.2, outlier.shape = NA) +
  ylim(c(0,70))+
  ylab('max speed (mile/hour)')+
  facet_grid(~geotype)+
  ggtitle('NONCBSA max speed from parabolic fit')
ggsave(file.path(figuredir, 'rural_ntype.maxspd.parabolic.by.cluster.pdf'), height =3, width = 7.5)

