# initialization of the clustering
library('pacman')

p_load(dplyr, MASS, reshape2,cowplot,rgdal, # rgdal: reading shapefiles
       ggplot2, corrplot, raster,readxl,RColorBrewer,fmsb) # sp: reading shapefiles; RColorBrewer color paletters; fmsb: radar charts
p_load(ncdf4,rgdal,ggmap,lawn,sp,shapefiles, maps, sf, fields, Imap, raster,readxl,
       tictoc,stargazer,psych,GPArotation,spdplyr,sp,shapefiles,tmap) # tictoc: time procedures
p_load(cluster,factoextra,DandEFA, xtable,psychTools, # DandEFA: for dandelion plots of EFA factors; factoextra: metrics to choose number of clusters
       xtable,psychTools,aCRM,clusterCrit,data.table,tigris,DAAG, fastDummies,tidyr,mlogit, filling, FNN)
       # clusterCrit: for cluster internal and external validity metrics; aCRM: for data imputation; FNN: knn regression
                 
