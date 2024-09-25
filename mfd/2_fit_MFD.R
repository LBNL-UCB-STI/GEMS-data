# 8/15/2024 LJ
# Rscript to fit the upperbound with a number of MFD parameter estimation methods
# 



mywd = "C:/Users/xiaodanxu/Documents/GitHub/GEMS-data/mfd"
#mywd <- "/Users/lingjin/Dropbox/Research/FHWA_GeoType/Work/Rscripts"
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
p_load(rlist) # for multiple plots
p_load(ggpubr) # for multiple plots

#source('initialization.R')
source('functions.R')
auxdatdir <- 'C:/data/CATTLab_delivery/AuxiliaryData' # the blob-microtype ids
rdatadir <-  "C:/data/RData/MFD_v2"
figuredir <- "C:/data/Figures/MFD_v2"
#outdir  <-  "C:/data/Data/"
part.dir <- 'C:/data/CATTLab_delivery/National/proc_data/partition_input' # partition xwalk


cat.downloads <- 'C:/data/CATTLab_delivery/National/'

cat.proc = 'C:/data/CATTLab_delivery/National/proc_data/all_tracts_frc3_5'  # tract level data aggregating frc3-5 consistently
cat.proc2 <- 'C:/data/CATTLab_delivery/National/proc_data/m3_tracts_frc1_5' # only tract level data for microtype 3 aggregating frc1-5

# 1. loop over states and fit and accumulate results by tracts
upfiles <- list.files(cat.proc, pattern = 'bytract_upperbounds_state')

fit.rural.all = NULL
fit.parabolic.all = NULL
fit.spddens.all = NULL

for(fnm in upfiles){
  dat = fread(file.path(cat.proc, fnm)) %>%
    filter(agg_density_plpm >= 1) %>% # filter out very low density and super-high speed
    group_by(tract_geoid)%>%
    mutate(nsample = n()) %>% # only esitmate when there are enough observations >= 100, 
    filter(nsample >=100)%>%
    select(-nsample)%>%
    ungroup()
  
  mfd.res_rural = dat %>%
    group_by(tract_geoid, microtype,state)%>%
    do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'rural'))
  
  mfd.res_urban = dat %>%
    group_by(tract_geoid, microtype,state)%>%
    do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'urban'))

  mfd.res_spddens = dat %>%
    group_by(tract_geoid, microtype,state)%>%
    do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'speed-density'))
  
  fit.rural.all = rbind(fit.rural.all, mfd.res_rural)
  fit.parabolic.all = rbind(fit.parabolic.all, mfd.res_urban)
  fit.spddens.all = rbind(fit.spddens.all, mfd.res_spddens)
  
  
}

save(fit.rural.all, fit.parabolic.all, fit.spddens.all, file = file.path(cat.proc, 'tract.mfd.fit.res.RData'))

# 
# 
# 
# mfd.func <- function(density){return(freeflow_speed*exp(-density/critical_density)*density)}
# ggplot(kk, aes(x = agg_density_plpm, y = agg_flow_plph )) +geom_point()+
#   geom_function(fun = mfd.func, geom = "line", color = 'red')


# 2. fit by partitions and get results

load( file = file.path(cat.proc, 'reaggreg_partition_upperbound.RData')) # res.upperbound,

dat = res.upperbound %>%
  filter(agg_density_plpm >= 1) %>% # filter out very low density and super-high speed
  group_by(partid)%>%
  mutate(nsample = n()) %>% # only esitmate when there are enough observations >= 100, 
  filter(nsample >=100)%>%
  select(-nsample)%>%
  ungroup()

mfd.res_steplinear = dat %>%
  group_by(partid, cbsa, cbsaname, state)%>%
  do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'rural'))

mfd.res_parabolic = dat %>%
  group_by(partid, cbsa, cbsaname, state)%>%
  do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'urban'))

mfd.res_spddens = dat %>%
  group_by(partid, cbsa, cbsaname, state)%>%
  do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'speed-density'))

save(mfd.res_steplinear, mfd.res_parabolic, mfd.res_spddens, file = file.path(cat.proc, "partition.mfd.fit.RData"))

############
# 3. fit the m3 frc1-5 by tract upper bounds
upfiles <- list.files(cat.proc2, pattern = 'upperbounds')

fitm3.rural.all = NULL
fitm3.parabolic.all = NULL
fitm3.spddens.all = NULL

for(fnm in upfiles){
  dat = fread(file.path(cat.proc2, fnm)) %>%
    filter(agg_density_plpm >= 1) %>% # filter out very low density and super-high speed
    group_by(tract_geoid)%>%
    mutate(nsample = n()) %>% # only esitmate when there are enough observations >= 100, 
    filter(nsample >=100)%>%
    select(-nsample)%>%
    ungroup()
  
  mfd.res_rural = dat %>%
    group_by(tract_geoid, microtype,state)%>%
    do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'rural'))
  
  mfd.res_urban = dat %>%
    group_by(tract_geoid, microtype,state)%>%
    do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'urban'))
  
  mfd.res_spddens = dat %>%
    group_by(tract_geoid, microtype,state)%>%
    do(fit.mfd.unified(input.density = .$agg_density_plpm ,input.flow = .$agg_flow_plph , mfd.shape = 'speed-density'))
  
  fitm3.rural.all = rbind(fitm3.rural.all, mfd.res_rural)
  fitm3.parabolic.all = rbind(fitm3.parabolic.all, mfd.res_urban)
  fitm3.spddens.all = rbind(fitm3.spddens.all, mfd.res_spddens)
  
  
}

save(fitm3.rural.all, fitm3.parabolic.all, fitm3.spddens.all, file = file.path(cat.proc2, 'tract.m3.mfd.fit.res.RData'))
# 
# 
# 
# # some diagnostics of results from alternative methods of fitting
# 
# ####################
# # merge the fit on partiions of 2+ tracts with other orphane/single tract fit in cbsas
# load(file = file.path(cat.proc, "partition.mfd.fit.RData")) # mfd.res_steplinear, mfd.res_parabolic, mfd.res_spddens, 
# 
# load(file = file.path(part.dir,'single_tract_partition_xwalk.RData')) # load pt1
# tmp.steplinear = fit.rural.all %>%
#   mutate(tract_geoid = as.numeric(tract_geoid)) %>%
#   filter(tract_geoid %in% pt1$tract_geoid)
# tmp.steplinear = tmp.steplinear%>%
#   left_join(pt1) %>%
#   ungroup() %>%
#   select(partid, cbsa, cbsaname,state, ff.spd.steplinear : R2.fit.steplinear)
# 
# tmp.parabolic = fit.parabolic.all %>%
#   mutate(tract_geoid = as.numeric(tract_geoid)) %>%
#   filter(tract_geoid %in% pt1$tract_geoid)%>%
#   left_join(pt1) %>%
#   ungroup() %>%
#   select(partid, cbsa, cbsaname,state, coef.x2.parabolic : R2.fit.parabolic)
# 
# tmp.spddens =  fit.spddens.all %>%
#   mutate(tract_geoid = as.numeric(tract_geoid)) %>%
#   filter(tract_geoid %in% pt1$tract_geoid)%>%
#   left_join(pt1) %>%
#   ungroup() %>%
#   select(partid, cbsa, cbsaname,state,freeflow_speed:R2.fit.speed_density)
# 
# tmp = full_join(tmp.steplinear, tmp.parabolic)
# tmp = full_join(tmp, tmp.spddens)  # mfd fit results all the orphan tract partitions
# 
# tmp$partitiontype <- 'orphan_tract'
# 
# ## also merge the fit on nonorphane partitions
# tmp2 <- full_join(mfd.res_steplinear, mfd.res_parabolic)
# tmp2<- full_join(tmp2, mfd.res_spddens)
# tmp2$partitiontype <- 'non_orphan_tract'
# 
# cbsa_fit_res <- bind_rows(tmp, tmp2)
# 
# save(cbsa_fit_res, file = file.path(rdatadir,'cbsa_partition_mfdfit.RData'))
# 
# # quantiles
# quantile(cbsa_fit_res$R2.fit.speed_density, seq(0.1, 1, 0.1))
# 
# ###############
# # merge with the tract_geoid then merge to the clustering results
# load(file = file.path(part.dir, 'all_partition_xwalk.RData')) # load pt
# load(file = file.path(part.dir, 'clustering_xwalk.RData')) # load all.clusters
# all.cluster = all.cluster %>%
#   select(is_cbsa, cbsa, cluster, tract_geoid20, tract_geoid)
# 
# tract.cbsa_fit_res <- cbsa_fit_res %>%
#   right_join(pt) %>%
#   rename(cbsa10 = cbsa,
#          cbsaname10 = cbsaname,
#          st_code10 = st_code)%>%
#   left_join(all.cluster)
# 
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster)), aes(x = cluster, y = R2.fit.speed_density, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   ylab('R2 from speed density') +
#   ggtitle('R2 from speed density')
# ggsave(file.path(figuredir, 'cbsa.r2.spdens.by.cluster.pdf'), height =3, width = 6.5)
#   
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster), max_capacity>0, max_capacity <2000), aes(x = cluster, y = max_capacity, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   ylim(c(0, 2500))+
#   ylab('max capacity from speed-density')+
#   ggtitle('max capacity from speed-density fit')
# ggsave(file.path(figuredir, 'cbsa.capacity.spdens.by.cluster.pdf'), height =3, width = 6.5)
# 
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster), max_capacity>0, max_capacity <2000), aes(x = cluster, y = critical_density, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2) +
#   ylim(c(0,300))+
#   ylab('critical_density from speed-density')+
#   ggtitle('critical_density from speed-density fit')
# ggsave(file.path(figuredir, 'cbsa.cridensity.spdens.by.cluster.pdf'), height =3, width = 6.5)
# 
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster), max_capacity>0, max_capacity <2000), aes(x = cluster, y = freeflow_speed, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   ylim(c(10,80))+
#   ylab('freeflow_speed from speed-density')+
#   ggtitle('freeflow_speed from speed-density fit')
# ggsave(file.path(figuredir, 'cbsa.ffspd.spdens.by.cluster.pdf'), height =3, width = 6.5)
# 
# 
# # plot the same for parabolic fit
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster)), aes(x = cluster, y = R2.fit.parabolic, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   ylab('R2 from parabolic') +
#   ggtitle('R2 from parabolic')
# ggsave(file.path(figuredir, 'cbsa.r2.parabolic.by.cluster.pdf'), height =3, width = 6.5)
# 
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster), max_capacity>0, max_capacity <2000), aes(x = cluster, y = capcity.fit.parabolic, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   ylim(c(0, 2500))+
#   ylab('max capacity from parabolic')+
#   ggtitle('max capacity from parabolic fit')
# ggsave(file.path(figuredir, 'cbsa.capacity.parabolic.by.cluster.pdf'), height =3, width = 6.5)
# 
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster), max_capacity>0, max_capacity <2000), aes(x = cluster, y = critical.density.fit.parabolic, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   ylim(c(0,300))+
#   ylab('critical_density from parabolic')+
#   ggtitle('critical_density from parabolic fit')
# ggsave(file.path(figuredir, 'cbsa.cridensity.parabolic.by.cluster.pdf'), height =3, width = 6.5)
# 
# ggplot(tract.cbsa_fit_res %>% filter(!is.na(cluster), max_capacity>0, max_capacity <2000), aes(x = cluster, y = max.spd.parabolic, color = cluster))+
#   geom_violin(width = 1) +
#   geom_boxplot(width = 0.2, outlier.shape = NA) +
#   ylim(c(10,80))+
#   ylab('freeflow_speed from parabolic')+
#   ggtitle('freeflow_speed from parabolic fit')
# ggsave(file.path(figuredir, 'cbsa.ffspd.parabolic.by.cluster.pdf'), height =3, width = 6.5)
# 
# 
# 
# ##################################################################
# # some plots #####################################################
# ##################################################################
# 
# # compare R2 and fitted parameters (capacity, critical density)
# ggplot(cbsa_fit_res, aes(x = R2.fit.parabolic, y = R2.fit.speed_density))+
#   geom_point()+
#   stat_density_2d(aes(fill = ..level..), geom = 'polygon') +
#   scale_fill_viridis_c(name = "density")+
#   geom_abline(intercept = 0, slope = 1, color = 'orange')+
#   xlim(c(0,1))+
#   ylim(c(0,1))+
#   theme_bw()+
#   ggtitle('R2 comparision')
# ggsave(file=file.path(figuredir,'R2.compare.cbsa.partitions.pdf'), height = 3, width = 4)
# 
# # capacity
# ggplot(cbsa_fit_res, aes(x = capcity.fit.parabolic, y = max_capacity))+
#   geom_point()+
#   stat_density_2d(aes(fill = ..level..), geom = 'polygon') +
#   scale_fill_viridis_c(name = "density")+
#   geom_abline(intercept = 0, slope = 1, color = 'orange')+
#   xlim(c(0,1500))+
#   ylim(c(0,1500))+
#   xlab('parabolic fit capacity')+
#   ylab('speed-desnity fit capacity')+
#   theme_bw()+
#   ggtitle('capacity comparision')
# ggsave(file=file.path(figuredir,'capacity.compare.cbsa.partitions.pdf'), height = 3, width = 4)
# 
# #with(cbsa_fit_res%>%filter(max_capacity<1000, capcity.fit.parabolic < 1000), cor(max_capacity, capcity.fit.parabolic ))
# # ~0.3
# 
# # critical density
# ggplot(cbsa_fit_res, aes(x = critical.density.fit.parabolic, y = critical_density))+
#   geom_point()+
#   stat_density_2d(aes(fill = ..level..), geom = 'polygon') +
#   scale_fill_viridis_c(name = "density")+
#   geom_abline(intercept = 0, slope = 1, color = 'orange')+
#   xlim(c(0,100))+
#   ylim(c(0,100))+
#   xlab('parabolic fit critical density')+
#   ylab('speed-desnity fit critical density')+
#   theme_bw()+
#   ggtitle('critical density comparision')  # however corr ~ 0.3
# ggsave(file=file.path(figuredir,'criticaldensity.compare.cbsa.partitions.pdf'), height = 3, width = 4)
# 
# #with(cbsa_fit_res%>%filter(critical_density<100, critical.density.fit.parabolic < 100), cor(critical.density.fit.parabolic, critical_density ))
# # ~0.2
# 
# 
# ###############3
# # select a few partitions to visualize the speed-density fit
# 
# set.seed(1395837)
# sel <- mfd.res_spddens %>%
#   ungroup() %>%
#   filter(R2.fit.speed_density>0.90) %>%   # 92 partitions with R2>0.95
#   sample_n(9)
# 
# # plot the 9 sample tracts and their fit
# figlist <- list()
# for(pid in sel$partid){
#   kk = res.upperbound %>%
#     filter(partid == pid) %>%
#     filter(agg_density_plpm  >=1)
#   freeflow_speed = sel$freeflow_speed[sel$partid == pid]
#   critical_density = sel$critical_density[sel$partid == pid]
#   r2.spd = round(sel$R2.fit.speed_density[sel$partid == pid], 2)
#   
#   a = mfd.res_parabolic$a.parabolic[mfd.res_parabolic$partid == pid]
#   b =  mfd.res_parabolic$b.parabolic[mfd.res_parabolic$partid == pid]
#   r2.para =  round(mfd.res_parabolic$R2.fit.parabolic[mfd.res_parabolic$partid == pid], 2)
#   # print(a)
#   # print(b)
#   
#   mfd.func <- function(density){return(freeflow_speed*exp(-density/critical_density)*density)}
#   mfd.func.par <- function(density){return(a*density*(density - b))}
#   
#   mfd.func.spd <- function(density){return(freeflow_speed*exp(-density/critical_density))}
#   
# # # for illustration
# #   ggplot(kk, aes(x = agg_density_plpm, y = agg_flow_plph )) +
# #     xlim(c(0,200))+
# #     ylim(c(0, 1000)) +
# #     geom_function(fun = mfd.func, geom = "line", color = 'red')
# #   ggsave(file = file.path(figuredir,'illustrate_flow-density.pdf'), height = 3, width = 4)
# # 
# #   ggplot(kk, aes(x = agg_density_plpm, y = aggspeed_mph )) +
# #     xlim(c(0,200))+
# #     ylim(c(0, 65)) +
# #     geom_function(fun = mfd.func.spd, geom = "line", color = 'red')
# #   ggsave(file = file.path(figuredir,'illustrate_speed-density.pdf'), height = 3, width = 4)
# #   
#   
#  p<- ggplot(kk, aes(x = agg_density_plpm, y = agg_flow_plph )) +geom_point()+
#     xlim(c(0,200))+
#     ylim(c(0, 2000)) +
#     geom_function(fun = mfd.func, geom = "line", color = 'red')+
#     geom_function(fun = mfd.func.par, geom = "line", color = 'blue')+
#    geom_vline(xintercept = critical_density, color = 'red', linetype = 2) +
#    geom_vline(xintercept = b/2, color = 'blue', linetype =2) +
#     theme_light()+
#    ggtitle(paste('partition id', pid, 'R2.red = ', r2.spd, "R2.blue =",r2.para ))
#  
#  p
#   ggsave(file = file.path(figuredir,paste0('r90.flowscatter_partition_',pid,'.pdf')), height = 3, width = 5)
#   
#  p2 <-  ggplot(kk, aes(x = agg_density_plpm, y = aggspeed_mph )) +geom_point()+
#    xlim(c(0,200))+
#    ylim(c(0, 65)) +
#    geom_function(fun = mfd.func.spd, geom = "line", color = 'red')+
#    geom_vline(xintercept = critical_density, color = 'red', linetype = 2) +
#    theme_light()+
#    ggtitle(paste('partition id', pid, 'R2.red = ', r2.spd ))
#  
#  ggarrange(p2, p, ncol = 2, nrow = 1)
#  ggsave(file = file.path(figuredir,paste0('r90.speedcatter_partition_',pid,'.pdf')), height = 3, width = 9)
#  
#   
# # figlist = list.append(figlist, p)  
# }
# 
# # ggarrange(plotlist = figlist,
# #           ncol = 3, nrow = 3)
# 
# 
# set.seed(1395837)
# sel <- mfd.res_spddens %>%
#   ungroup() %>%
#   filter(R2.fit.speed_density>0.80 & R2.fit.speed_density <0.90) %>%   
#   sample_n(9)
# 
# # plot the 9 sample tracts and their fit
# figlist <- list()
# for(pid in sel$partid){
#   kk = res.upperbound %>%
#     filter(partid == pid) %>%
#     filter(agg_density_plpm  >=1)
#   freeflow_speed = sel$freeflow_speed[sel$partid == pid]
#   critical_density = sel$critical_density[sel$partid == pid]
#   r2.spd = round(sel$R2.fit.speed_density[sel$partid == pid], 2)
#   
#   a = mfd.res_parabolic$a.parabolic[mfd.res_parabolic$partid == pid]
#   b =  mfd.res_parabolic$b.parabolic[mfd.res_parabolic$partid == pid]
#   r2.para =  round(mfd.res_parabolic$R2.fit.parabolic[mfd.res_parabolic$partid == pid], 2)
#   # print(a)
#   # print(b)
#   
#   mfd.func <- function(density){return(freeflow_speed*exp(-density/critical_density)*density)}
#   mfd.func.par <- function(density){return(a*density*(density - b))}
#   
#   mfd.func.spd <- function(density){return(freeflow_speed*exp(-density/critical_density))}
#   
#   
#   p<- ggplot(kk, aes(x = agg_density_plpm, y = agg_flow_plph )) +geom_point()+
#     xlim(c(0,200))+
#     ylim(c(0, 2000)) +
#     geom_function(fun = mfd.func, geom = "line", color = 'red')+
#     geom_function(fun = mfd.func.par, geom = "line", color = 'blue')+
#     geom_vline(xintercept = critical_density, color = 'red', linetype = 2) +
#     geom_vline(xintercept = b/2, color = 'blue', linetype =2) +
#     theme_light()+
#     ggtitle(paste('partition id', pid, 'R2.red = ', r2.spd, "R2.blue =",r2.para ))
#   
#   p
#   ggsave(file = file.path(figuredir,paste0('r80.flowscatter_partition_',pid,'.pdf')), height = 3, width = 5)
#   
#   p2 <-  ggplot(kk, aes(x = agg_density_plpm, y = aggspeed_mph )) +geom_point()+
#     xlim(c(0,200))+
#     ylim(c(0, 65)) +
#     geom_function(fun = mfd.func.spd, geom = "line", color = 'red')+
#     geom_vline(xintercept = critical_density, color = 'red', linetype = 2) +
#     theme_light()+
#     ggtitle(paste('partition id', pid, 'R2.red = ', r2.spd ))
#   
#   ggarrange(p2, p, ncol = 2, nrow = 1)
#   ggsave(file = file.path(figuredir,paste0('r80.speedcatter_partition_',pid,'.pdf')), height = 3, width = 9)
#   
#   
#   # figlist = list.append(figlist, p)  
# }
# 
# set.seed(1395837)
# sel <- mfd.res_spddens %>%
#   ungroup() %>%
#   filter(R2.fit.speed_density>0.70 & R2.fit.speed_density <0.80) %>%   
#   sample_n(9)
# 
# # plot the 9 sample tracts and their fit
# figlist <- list()
# for(pid in sel$partid){
#   kk = res.upperbound %>%
#     filter(partid == pid) %>%
#     filter(agg_density_plpm  >=1)
#   freeflow_speed = sel$freeflow_speed[sel$partid == pid]
#   critical_density = sel$critical_density[sel$partid == pid]
#   r2.spd = round(sel$R2.fit.speed_density[sel$partid == pid], 2)
#   
#   a = mfd.res_parabolic$a.parabolic[mfd.res_parabolic$partid == pid]
#   b =  mfd.res_parabolic$b.parabolic[mfd.res_parabolic$partid == pid]
#   r2.para =  round(mfd.res_parabolic$R2.fit.parabolic[mfd.res_parabolic$partid == pid], 2)
#   # print(a)
#   # print(b)
#   
#   mfd.func <- function(density){return(freeflow_speed*exp(-density/critical_density)*density)}
#   mfd.func.par <- function(density){return(a*density*(density - b))}
#   
#   mfd.func.spd <- function(density){return(freeflow_speed*exp(-density/critical_density))}
#   
#   
#   p<- ggplot(kk, aes(x = agg_density_plpm, y = agg_flow_plph )) +geom_point()+
#     xlim(c(0,200))+
#     ylim(c(0, 2000)) +
#     geom_function(fun = mfd.func, geom = "line", color = 'red')+
#     geom_function(fun = mfd.func.par, geom = "line", color = 'blue')+
#     geom_vline(xintercept = critical_density, color = 'red', linetype = 2) +
#     geom_vline(xintercept = b/2, color = 'blue', linetype =2) +
#     theme_light()+
#     ggtitle(paste('partition id', pid, 'R2.red = ', r2.spd, "R2.blue =",r2.para ))
#   
#   p
#   ggsave(file = file.path(figuredir,paste0('r70.flowscatter_partition_',pid,'.pdf')), height = 3, width = 5)
#   
#   p2 <-  ggplot(kk, aes(x = agg_density_plpm, y = aggspeed_mph )) +geom_point()+
#     xlim(c(0,200))+
#     ylim(c(0, 65)) +
#     geom_function(fun = mfd.func.spd, geom = "line", color = 'red')+
#     geom_vline(xintercept = critical_density, color = 'red', linetype = 2) +
#     theme_light()+
#     ggtitle(paste('partition id', pid, 'R2.red = ', r2.spd ))
#   
#   ggarrange(p2, p, ncol = 2, nrow = 1)
#   ggsave(file = file.path(figuredir,paste0('r70.speedcatter_partition_',pid,'.pdf')), height = 3, width = 9)
#   
#   
#   # figlist = list.append(figlist, p)  
# }
# 
# 
# 
# 
