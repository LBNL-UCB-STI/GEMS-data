# Master File to generate microtypes and geotypes 
# from the cleaned input data
# NATALIE POPOVICH
# BERKELEY NATIONAL LAB
# MAIN FILE: SEP 16 2020
# UPDATED: Jan 12 2021
####################################3

# set working directory and sub-directories
mywd <- "~/Library/CloudStorage/Box-Box/FHWA/Task1/TransportGeography-Revision"
mywd <- '/home/ubuntu/Task1/TransportGeography-Revision' # if on AWS
mywd <- 'C:/FHWA/For FHWA folks/TransportGeography-Revision' # if on data bucket

setwd(mywd)

figuredir <- "./Figures"
inputsdir <- "./Data/InputData"
datadir <- "./Data/OutputData"
rdatadir <-  "./Data/RData"
tabdir <- './Tables'
resultsdir <- './Results'

# Load prep files
source('initialization.R')
source('functions.R')
################
# Code to generate microtypes, including figures
###############
source('1_transgeoV2_stage1_prep.R') # checked
source('2_transgeoV2_stage1_clustering_v4.R') # checked
    # Figure 2: factor loadings 
    # Figure 3: spider plots of microtypes

########
# robustness checks for the microtypes
###########
# microtypes
source('3_transgeoV2_stage1_clusterValidation.R') #checked

################
# Code to generate geotypes, including figures and tables
###############
source('4_transgeoV2_stage2_compile_v2.R') #DONE
source('5_transgeoV2_stage2_clustering_v7.R') #checked
source('6_transgeoV2_stage2_plots.R') #checked

############# 
# Mapping results as layers for Mapbox
############
#source('7_transgeoV2_maps.R') # DONE


# Additional code not in the publication, but used for the GEMS model
###############
# source()


