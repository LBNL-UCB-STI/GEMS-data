###############
## CODE FOR FIRST-STAGE CLUSTERING TO GENERATE MICROTYPES
        # IMPORT CLEANED DATASET OF RAW INPUTS
        # GENERATE FACTORS WITH EXPLORATORY FACTOR ANALYSIS
        # CLUSTER USING CLARA
        # SAVE SHAPEFILES OF MICROTYPES

# NATALIE POPOVICH
# LAWRENCE BERKELEY NATIONAL LAB
# LAST UPDATED: DEC 19 2019
# modified Ling Jin: 2/24/2020
# adding variations and compare results
################
# this version removes all county-level variables (weather and macro-econ indicators)
# set working directory
#mywd <- "/Users/lingjin/Dropbox/Research/FHWA_GeoType/Work/Rscripts" #Ling2
# #mywd = "C:/Users/ljin/Dropbox/Research/FHWA_GeoType/Work/Rscripts" #Ling 1
# mywd <- "~/Box/FHWA/" # Natalie
# setwd(mywd)
# 
# figuredir <- "./Figures/TransGeo"
# #datadir1  <-  "../Data/JGEOTRANS_PaperRevision" # Ling's CPU
# datadir <- "./Data/CleanData/" # houses most of the data non-specific to TransGEo paper
# datadir1 <- "./Data/CleanData/TransGeo" # houses data specific to TransGeo methods
# rdatadir <-  "./RData/TransGeo"
# tabdir <- './Tables/TransGeo'
# resultsdir <- './Results/TransGeo'
# 
# source('initialization.R')
# source('functions.R')

####################################################
# prepare for clustering comparison

factors <- c("MR1", "MR2", "MR3","MR4", "MR5", "MR6", "MR7", "MR8", "MR9", "MR10","MR11","MR12" )

# load dataset of factors and loadings. 
# clustering datasets
load(file = file.path(rdatadir,'factor12.scores.RData')) # load scores of the 12 factors: cluster.df
load(file =  file.path(rdatadir,'raw.scaled.RData')) # load the scaled raw data: raw.df

####################################################
# CLUSTERING robustness check
####################################################

# cluster quality using raw data vs data of reduced dimension. The reduced dimension data is supposed to produce better quality clusters because we removed the random variations and mitigate the "curse of dimensionality"
# as it's memory intensive, use bootstrap to exhaust data space
 
# geotypes are smaller size data, should increase the draw, try draw in 1:100 to get better accuracy of the final esitmates. 

# cluster quality metric: Silhouette and dbi

# some burden tests showed that samples = 1000 and sampsize 100 (note sampsize 100 is greater than the default) 
# can improve the speed while produce similar results.
set.seed(1) # so results are replicable each time

tic()
cluster6 <- clara(cluster.df, 6, metric = "euclidean", #euclidean distance metric
                  stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
toc()

set.seed(1) # so results are replicable each time
tic()
test <- clara(cluster.df, 6, metric = "euclidean", #euclidean distance metric
              stand = FALSE, samples = 2000, sampsize = 70, pamLike = TRUE) #partition around mediods
toc()

table(test$clustering, cluster6$clustering)

#gg = fviz_nbclust(sdat, cluster::clara,  method = "silhouette",  k.max = 10) 
#plot(gg)

tic()
reduced.res = NULL

for(draw in 1:20){ # bootstrap 20 times
  
  set.seed(draw)
  sdat = sample_n(cluster.df, 30000) # reduced data
#  print(sdat[1,])
  print(paste('draw =',draw))
  
  for(k in 2:10){
    print(paste('k=',k))
    # using reduced dimension data
    cluster.sdat <- clara(sdat, k, metric = "euclidean", #euclidean distance metric
                     stand = FALSE, samples = 1000, sampsize = 100, pamLike = TRUE) #partition around mediods
    
    indexVal.sdat = intCriteria(as.matrix(sdat),
                      cluster.sdat$clustering,c("Davies_Bouldin","Silhouette"))
     
    tmp = data.frame(k = rep(k,2),
                     index.value = c(indexVal.sdat$davies_bouldin, 
                                     indexVal.sdat$silhouette),
                     draw.n = rep(draw,2),
                     index.name = c('davies_bouldin','silhouette'))
    
    reduced.res = rbind(reduced.res,tmp)

    rm(tmp);rm(indexVal.sdat)
    
  }
}
  
toc()

save(reduced.res, file = file.path(rdatadir,'dbi.asw.fact12data.stage1.RData'))

reduced.res[reduced.res$index.name=='davies_bouldin','index.value'] = 
  1/reduced.res[reduced.res$index.name=='davies_bouldin','index.value']
# remove outliers
reduced.res = reduced.res %>%
  filter(index.value <3)



ggplot(reduced.res, 
       aes(x = as.factor(k), y = index.value, fill = index.name)) +
  geom_boxplot()


dat = reduced.res  %>%
  group_by(k,index.name) %>%
  summarize(mean = mean(index.value),
            sd = sd(index.value),
            median = median(index.value),
            q25 = quantile(index.value,0.25),
            q75 = quantile(index.value,0.75))

ggplot(dat,aes(x= as.integer(k),y=mean,color=index.name))+
  geom_rect(aes(xmin=4, xmax=6, ymin=-Inf, ymax=Inf),
            fill = 'lightgrey',alpha = 0.8,inherit.aes = F, color = 'white')+
  geom_point()+
    geom_line() +
    geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.3) +
    xlab('Number of Clusters') +
  ylab('Cluster Validity Metrics')+
    scale_x_continuous(breaks=2:10)+
    scale_y_continuous(limits = c(0,0.8)) +
    scale_color_discrete(name = "Validity Metrics", labels = c("DBI", "ASW"))  +
  theme_bw() 


 ggsave(file = file.path(figuredir,'stage1_cluster_validity.pdf'),
        height = 4, width = 6) 
 
 dat1 = dat %>%
   filter(index.name == 'davies_bouldin' )
 
 dat2 = dat %>%
   filter(index.name == 'silhouette')
 


 ggplot(dat1, aes(x = k, y= mean)) +
   theme_bw()+
   geom_vline(xintercept = 6, size = rel(1.1), color = 'grey') +
   geom_point(color = '#f59ff2') +
   geom_line(color = '#f59ff2') +
   geom_point(data = dat2, aes(x = k, y = mean*4), color = "#7be2ed")+
   geom_line(data = dat2, aes(x = k, y = mean*4),color = "#7be2ed") +
   scale_y_continuous(
     # Features of the first axis
     name = "Inverse DBI",
     # Add a second axis and specify its features
     sec.axis = sec_axis(~./4, name="Silhouette")
   ) + 

   theme(
     axis.title.y = element_text(color ='#f59ff2', size=13),
     axis.title.y.right = element_text(color = "#7be2ed", size=13)) +
   xlab('Number of Clusters') 
 
 ggsave(file = file.path(figuredir,'stage1_cluster_validity_pick6.pdf'),
        height = 3, width = 5) 

 
 ggplot(dat,aes(x= as.integer(k),y=median,color=index.name))+
   theme_bw() +
   geom_rect(aes(xmin=4, xmax=6, ymin=-Inf, ymax=Inf),
             fill = 'lightgrey',alpha = 0.8,inherit.aes = F, color = 'white')+
   geom_point()+
   geom_line() +
   geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.3) +
   xlab('Number of Clusters') +
   ylab('Cluster Validity Metrics')+
   scale_x_continuous(breaks=2:10)+
   scale_y_continuous(limits = c(0,0.8)) +
   scale_color_discrete(name = "Validity Metrics", labels = c("DBI", "ASW"))
 

# using raw data
tic()
raw.res = NULL
for(draw in 1:20){ # bootstrap 20 times
  
  set.seed(draw)
  rdat = sample_n(raw.df,30000)
  print(paste('draw =',draw))
  
  for(k in 2:10){
    print(paste('k=',k))
    
    cluster.rdat <- clara(rdat, k, metric = "euclidean", #euclidean distance metric
                          stand = FALSE, samples = 1000, sampsize = 100, pamLike = TRUE) #partition around mediods
    
    indexVal.rdat = intCriteria(as.matrix(rdat),
                                cluster.rdat$clustering,c("Davies_Bouldin","Silhouette"))
    
    tmp = data.frame(k = rep(k,2),
                     index.value = c(indexVal.rdat$davies_bouldin, 
                                     indexVal.rdat$silhouette),
                     draw.n = rep(draw,2),
                     index.name = c('davies_bouldin','silhouette'))
    
    raw.res = rbind(raw.res,tmp)
    
    rm(tmp);rm(indexVal.rdat)
    
  }
}

toc()

save(raw.res, file = file.path(rdatadir,'dbi.asw.rawdata.stage1.RData'))

asw.res = NULL
raw.res[raw.res$index.name=='davies_bouldin','index.value'] = 
  1/raw.res[raw.res$index.name=='davies_bouldin','index.value']
# remove outliers
raw.res = raw.res %>%
  filter(index.value <3)



ggplot(raw.res, 
       aes(x = as.factor(k), y = index.value, fill = index.name)) +
  geom_boxplot()




  
  
  
  
  gg = fviz_nbclust(sdat, cluster::clara,  method = "silhouette",  k.max = 10) 
  
  tmp = gg$data
  print(tmp)
  tmp$draw.n = draw
  tmp$index = 'asw'
  colnames(tmp) <- c('k','index.value','draw.n','index.name') 
  asw.res = rbind(asw.res,tmp)
  
  rm(tmp)
  rm(sdat)

toc()
save(asw.res,file = file.path(rdatadir,'asw.boot.factor12.clusters.RData'))


tic()
asw.res.raw = NULL
for(draw in 1:30){
  set.seed(draw)
  sdat = sample_n(raw.df, 30000)
  gg = fviz_nbclust(sdat, cluster::clara,  method = "silhouette",  k.max = 10) 
  
  tmp = gg$data
  tmp$draw.n = draw
  tmp$index = 'asw'
  colnames(tmp) <- c('k','index.value','draw.n','index.name') 
  asw.res.raw = rbind(asw.res.raw,tmp)
  
  rm(tmp)
  rm(sdat)
}
toc()

save(asw.res.raw,file = file.path(rdatadir,'asw.boot.rawdat.clusters.RData'))



tic()
dbi.res = NULL
for(draw in 1:20){
  set.seed(draw)
  sdat = sample_n(cluster.df, 30000)
  print(sdat[1,])
  print(draw)
  
  res = data.frame(k = 2:10, index.value = 2:10, index= rep('dbi',9))
  for(k in 2:10){
    cluster <- clara(sdat, k, metric = "euclidean", #euclidean distance metric
                     stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
    
    res[res$k==k,'index.value'] =intCriteria(as.matrix(sdat),cluster$clustering,"Davies_Bouldin")
  }
  res$draw.n = draw
  
  print(res)
  
  dbi.res = rbind(dbi.res,res)
  
  rm(res)
  rm(sdat)
}
toc()

save(dbi.res,file = file.path(rdatadir,'dbi.boot.factor12.clusters.RData'))


tic()
dbi.res.raw = NULL
for(draw in 1:20){
  set.seed(draw)
  sdat = sample_n(raw.df, 30000)
  print(sdat[1,])
  print(draw)
  
  res = data.frame(k = 2:10, index.value = 2:10, index= rep('dbi',9))
  for(k in 2:10){
    cluster <- clara(sdat, k, metric = "euclidean", #euclidean distance metric
                     stand = FALSE, samples = 1000,sampsize = 100, pamLike = TRUE) #partition around mediods
    
    res[res$k==k,'index.value'] =intCriteria(as.matrix(sdat),cluster$clustering,"Davies_Bouldin")
  
    
    }
  
  

  print(res)
  
  dbi.res.raw = rbind(dbi.res.raw,res)
  
  rm(res)
  rm(sdat)
}
toc()

save(dbi.res.raw,file = file.path(rdatadir,'dbi.boot.rawdat.clusters.RData'))



















##### computation burden tests###################
tic()
res = data.frame(k = 2:10, asw = 2:10)
for(k in 2:10){
  cluster <- clara(s1, k, metric = "euclidean", #euclidean distance metric
                    stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
  
  res[res$k==k,'asw'] =intCriteria(as.matrix(s1),cluster$clustering,"Davies_Bouldin")
}
toc()

ggplot(res,aes(x=k,y=asw))+geom_line()


tic()
cluster6 <- clara(s1, 6, metric = "euclidean", #euclidean distance metric
                  stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
toc() # 15 s

tic()
cluster6 <- clara(cluster.df, 6, metric = "euclidean", #euclidean distance metric
                  stand = FALSE, samples = 5000, pamLike = TRUE) #partition around mediods
toc() # 33 s


tic()
db =intCriteria(as.matrix(cluster.df),cluster6$clustering,"Davies_Bouldin")
toc() #~35 s for 30k sample, and 210 s or 3.6 min to run for full sample, the results are similar 1.327 vs 1.335

tic()
pbc = intCriteria(as.matrix(cluster.df),cluster6$clustering,"Point_Biserial")
toc() #~35 s for 30k sample, and 210 s or 3.6 min to run for full sample, the results are similar -0.64 vs Na, both are not reasonable

tic()
asw = intCriteria(as.matrix(cluster.df),cluster6$clustering,"Silhouette")
toc() #~35 s for 30k sample, and 210 s or 3.6 min to run for full sample, the results are similar 0.136 vs 0.126


