library(sf)
library(dplyr)
library(data.table)
library(tidycensus)
library(tools)

readRenviron("~/.Renviron")

merge_census_hpms<- function(state='california2017',state_tracts,year){
  print(state)
  hpms_geometry <- st_read(paste0(state,'/',toTitleCase(substr(state,10,50)), '.shp')) # read HPMS data
  crs_hpms <- st_crs(hpms_geometry)
  
  state_tracts <- st_transform(state_tracts, crs = crs_hpms)
  sf::sf_use_s2(FALSE) 
  print(nrow(hpms_geometry))
  # assign tracts to each link
  hpms_geometry_by_tracts = st_intersection(st_zm(hpms_geometry), state_tracts)
  print(nrow(hpms_geometry_by_tracts))
  hpms_geometry_by_tracts$Length = st_length(hpms_geometry_by_tracts) # re-generate link length after split network by tracts
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% mutate(lanemiles = as.numeric(Through_La * Length / 1609.34)) # compute lane miles
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% filter(lanemiles > 0) # remove invalid links
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% mutate(g_type = st_geometry_type(.)) 
  hpms_geometry_by_tracts <- hpms_geometry_by_tracts %>% filter(g_type %in% c('LINESTRING', 'MULTILINESTRING'))

  st_write(hpms_geometry_by_tracts, paste0( state, '_HPMS_with_',year,'_GEOID_LANEMILE.geojson'),append=FALSE, na='NA') # save merged data
}

# load HPMS data
setwd("C:/FHWA_R2/Network/RawData/") # set input directory
file_list <- list.files("HPMS2017", full.names = TRUE) # list all HPMS data files
file_list <- file_list[!grepl("zip", file_list, ignore.case = TRUE)]
file_list <- file_list[!grepl("test", file_list, ignore.case = TRUE)]
file_list <- file_list[!grepl("GEOID", file_list, ignore.case = TRUE)]

# load state census tracts (whole U.S.)
tracts <- st_read("C:/FHWA_R2/spatial_boundary/CleanData/combined_tracts_2020.geojson")

for (f in file_list){
  merge_census_hpms(f,tracts,2020)
}


# Produce CA tracts and hpms data for mapping
ca_tracts <- tracts %>% filter(STATEFP=='06')
st_write(ca_tracts, 'ca_tracts.geojson',append=FALSE, na='NA') # save filtered data

hpms_df = read.csv("C:/FHWA_R2/Network/network_microtype_metrics.csv")

ca_hpms_df <- hpms_df%>% filter(state=="california")
ca_hpms_df_slim <- ca_hpms_df[,-3]
ca_hpms_df_slim <- ca_hpms_df_slim[,-3]
# Convert both columns to character and add leading zeros
ca_tracts$geoid <-sprintf("%010s", ca_tracts$GEOID)
ca_hpms_df_slim$geoid <- paste0("0",ca_hpms_df_slim$tract) #sprintf("%010i", ca_hpms_df$tract)

merged_data <- merge(ca_tracts, ca_hpms_df_slim, by.x = "geoid", by.y = "geoid")

st_write(merged_data[,-1], 'ca_tracts_hpms_data4.geojson',append=FALSE, na='NA') # save merged data


# Produce lanemiles by F_System table
file = "HPMS2017/california2017"
hpms_geometry17 <- st_read(paste0(file,'/',toTitleCase(substr(file,10,50)), '.shp'))
hpms_geometry17$Length = st_length(hpms_geometry17) # re-generate link length after split network by tracts
hpms_geometry17 <- hpms_geometry17 %>% mutate(lanemiles = as.numeric(Through_La * Length / 1609.34)) # compute lane miles
#hpms_geometry17 <- hpms_geometry17 %>% filter(lanemiles > 0)
#print(aggregate(hpms_geometry17$lanemiles, by=list(Category= hpms_geometry17$F_System), FUN=sum))
print(nrow(hpms_geometry17))
print(sum(is.na(hpms_geometry17$F_System)))
print(sum(is.na(hpms_geometry17$IRI)))
print(sum(is.na(hpms_geometry17$Access_Con)))
print(sum(is.na(hpms_geometry17$Speed_Limi)))
print(sum(is.na(hpms_geometry17$AADT)))
print(sum(is.na(hpms_geometry17$Through_La)))
print(nrow(subset(hpms_geometry17,F_System>0)))
print(nrow(subset(hpms_geometry17,IRI>0)))
print(nrow(subset(hpms_geometry17,Access_Con>0)))
print(nrow(subset(hpms_geometry17,Speed_Limi>0)))
print(nrow(subset(hpms_geometry17,AADT>0)))
print(nrow(subset(hpms_geometry17,Through_La>0)))

print(nrow(filter(hpms_geometry17, (F_System>0&F_System<4)|NHS>0)))
print(nrow(filter(hpms_geometry17, ((F_System>0&F_System<4)|NHS>0)&IRI>0)))
print(nrow(filter(hpms_geometry17,  ((F_System>0&F_System<4)|NHS>0)&IRI==0)))
print(sum(is.na(filter(hpms_geometry17,  (F_System>0&F_System<4)|NHS>0)$IRI)))
print(nrow(filter(hpms_geometry17, F_System>1,F_System<4)))
print(nrow(filter(hpms_geometry17, F_System>1,F_System<4,Access_Con>0)))
print(nrow(filter(hpms_geometry17, F_System>1,F_System<4,Access_Con==0)))
print(sum(is.na(filter(hpms_geometry17, F_System>1 & F_System<4)$Access_Con)))

file = "HPMS2019/california2019"
hpms_geometry <- st_read('HPMS2019/california2019/california_pr_2019.shp')
hpms_geometry$Length = st_length(hpms_geometry) # re-generate link length after split network by tracts
hpms_geometry <- hpms_geometry %>% mutate(lanemiles = as.numeric(through_la * Length / 1609.34)) # compute lane miles
hpms_geometry <- hpms_geometry %>% filter(lanemiles > 0)
print(aggregate(hpms_geometry$lanemiles, by=list(Category= hpms_geometry$f_system), FUN=sum))


file = "HPMS2022/california2022"
hpms_geometry <- st_read(paste0(file,'/',toTitleCase(substr(file,10,50)), '.shp'))
hpms_geometry$Length = st_length(hpms_geometry) # re-generate link length after split network by tracts
hpms_geometry <- hpms_geometry %>% mutate(lanemiles = as.numeric(through_la * Length / 1609.34)) # compute lane miles
hpms_geometry <- hpms_geometry %>% filter(lanemiles > 0)
print(aggregate(hpms_geometry$lanemiles, by=list(Category=hpms_geometry$f_system), FUN=sum))
print(nrow(hpms_geometry))
print(sum(is.na(hpms_geometry$f_system)))
print(sum(is.na(hpms_geometry$iri)))
print(sum(is.na(hpms_geometry$access_con)))
print(sum(is.na(hpms_geometry$speed_limi)))
print(sum(is.na(hpms_geometry$aadt)))
print(sum(is.na(hpms_geometry$through_la)))
print(nrow(subset(hpms_geometry,f_system>0)))
print(nrow(subset(hpms_geometry,iri>0)))
print(nrow(subset(hpms_geometry,access_con>0)))
print(nrow(subset(hpms_geometry,speed_limi>0)))
print(nrow(subset(hpms_geometry,aadt>0)))
print(nrow(subset(hpms_geometry,through_la>0)))
print(nrow(subset(hpms_geometry,f_system==0)))
print(nrow(filter(hpms_geometry, ((f_system>0&f_system<4)| nhs>0))))
print(nrow(filter(hpms_geometry, ((f_system>0&f_system<4)| nhs>0) &iri>0)))
print(nrow(filter(hpms_geometry, ((f_system>0&f_system<4)| nhs>0) & iri==0)))
print(sum(is.na(filter(hpms_geometry, (f_system>0&f_system<4)| nhs>0)$iri)))
print(nrow(filter(hpms_geometry, f_system>1,f_system<4)))
print(nrow(filter(hpms_geometry, f_system>1,f_system<4,access_con>0)))
print(nrow(filter(hpms_geometry, f_system>1,f_system<4,access_con==0)))
print(sum(is.na(filter(hpms_geometry, f_system>1 & f_system<4)$access_con)))
write.csv(hpms_geometry, "HPMS2022/california2022/california_hpms2022.csv", row.names=FALSE)
