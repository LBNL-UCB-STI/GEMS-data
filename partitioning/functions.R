# global variables and functions used
# Ling Jin 
# last updated 5/19/2023

stage1.factor <- function(data_scaled, nfactor,cutoff = 0.3){
  # factor analysis in stage 1
  
  factor <- fa(data_scaled, nfactors = nfactor, # number of factors
                 rotate = "oblimin", # how (if at all) to rotate the factors
                 fm="minres", # what factor extraction method to use. here maximum likelihood
                 max.iter = 10000, missing = TRUE, impute = "median", warnings = TRUE)
#  print(factor)
#  print(factor$loadings,cutoff = cutoff) # show factor loadings and variation explained
  
#  loadings <- factor$loadings
  fa2latex(factor, digits=3,rowlabels=TRUE,apa=TRUE,short.names=FALSE,cumvar=T, cut=cutoff,big=.3,alpha=.05,font.size ="scriptsize",
           caption= paste('TransGeo',nfactor,"factor solution"),label="default",silent=FALSE,
           file= file.path(tabdir,paste0("transgeofac_loadings",nfactor,".tex")))
  return(factor)
  
}

names_ordered <- c("Pct water","Development Intensity",
                   "Avg.Circuity", "Dead-end Proportion", "Intersection Density",
                   "Self-loop Proportion", "Street Density", "Avg Street Length",
                   "Broadband", 
                   "Jobs-Housing Balance",
                   "Road grade", "Avg IRI",
                   "Non-attainment pollutants", 
                   "Pct Full Access Control","Pct Partial Access Control", 
                   "Pct local roads", "Pct Midsize Roads",
                   "Pct Highways" ,"Pct Truck AADT",
                   "Truck AADT per Lane Mile",
                   "Lane-miles per Sq. Km" ,
                   "Lane-meters per Capita", 
                   "Population Density", 
                   "Job Density",
                   "Pct manufacturing jobs", 
                   "Pct mining jobs",  
                   "Land use - Ag" )

transgeo_names_ordered = c("Broadband", "Non-attainment pollutants", "Land use - Ag", "Pct water",
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

transgeo_names_ordered2 = c("Broadband", "Non-attainment pollutants", "Land use - Ag", "Pct water",
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
                           'Trip Source Magnitude',
                           'Demand Supply Burden')


plot.loadings <- function(fa_output, factor.names = NULL, names_ordered = names_ordered){
  # fa_output is the output from fa()
  # factor.names is a character vector that user can supply to rename the factors
  # use 0.25 as cutoff
  require(fields)
  loadings <- fa_output$loadings
  dat = as.matrix(data.frame(matrix(as.numeric(loadings), attributes(loadings)$dim, dimnames=attributes(loadings)$dimnames)))
  # name the factors
  if(!is.null(factor.names)){
    colnames(dat) <- factor.names
  }
  
  mycols = colorRampPalette(c(tim.colors(10)[c(1,2,3)],'white','white',tim.colors(10)[c(8,9,10)]))

  tick<- seq(-1,1,0.25)
  par(mar=c(5,10,8,2.8) + 0.1)

  plot.dat = dat[rev(rownames(dat)),]

  image(1:dim(plot.dat)[2],1:dim(plot.dat)[1],round(t(as.matrix(plot.dat)),2),col = rev(mycols(8)),zlim=c(-1,1),breaks = seq(-1,1,0.25),
      yaxt='n',xaxt='n',
      xlab='',ylab='',cex.lab=1.2,cex.axis=1.4)
  abline(v=1:(dim(plot.dat)[2])+0.5,lwd=1.5)
  box(lwd=1.5)
  #axis(3, labels = FALSE)
  ## Create some text labels
  labels <- colnames(dat)
  ## Plot x axis labels at default tick marks
  text(x=1:(dim(plot.dat)[2])-0.7,y=rep((dim(plot.dat)[1]+1),dim(plot.dat)[2]), adj = 0, pos = 4,offset = 1, srt = 50, 
       labels = labels, xpd = TRUE,cex=0.8)
  
  #mtext(colnames(dat),3,at=1:10,las=2,adj=0,cex=0.8,line = 0.5,srt = 45)
  mtext(rev(names_ordered),2,at=1:(dim(plot.dat)[1]),las=2,adj=1,cex=0.85,line = 0.5)
  abline(h=1:(dim(plot.dat)[1])+0.5,lwd=0.3)
  
}

plot.loadings2 <- function(fa_output, factor.names = NULL, names_ordered = names_ordered){
  # fa_output is the output from fa()
  # factor.names is a character vector that user can supply to rename the factors
  # use 0.3 as cutoff
  require(fields)
  loadings <- fa_output$loadings
  dat = as.matrix(data.frame(matrix(as.numeric(loadings), attributes(loadings)$dim, dimnames=attributes(loadings)$dimnames)))
  # name the factors
  if(!is.null(factor.names)){
    colnames(dat) <- factor.names
  }
  
  mycols = colorRampPalette(c(tim.colors(10)[c(1,2)],'white','white',tim.colors(10)[c(9,10)]))
  
  tick<- c(-1,-0.65,-0.3,0,0.3,0.65,1)
  par(mar=c(5,10,8,2.8) + 0.1)
  
  plot.dat = dat[rev(rownames(dat)),]
  
  image(1:dim(plot.dat)[2],1:dim(plot.dat)[1],
        round(t(as.matrix(plot.dat)),2),col = rev(mycols(6)),zlim=c(-1,1),
        breaks = tick,
        yaxt='n',xaxt='n',
        xlab='',ylab='',cex.lab=1.2,cex.axis=1.4)
  abline(v=1:(dim(plot.dat)[2])+0.5,lwd=1.5)
  box(lwd=1.5)
  #axis(3, labels = FALSE)
  ## Create some text labels
  labels <- colnames(dat)
  ## Plot x axis labels at default tick marks
  text(x=1:(dim(plot.dat)[2])-0.7,y=rep((dim(plot.dat)[1]+1),dim(plot.dat)[2]), adj = 0, pos = 4,offset = 1, srt = 50, 
       labels = labels, xpd = TRUE,cex=0.8)
  
  #mtext(colnames(dat),3,at=1:10,las=2,adj=0,cex=0.8,line = 0.5,srt = 45)
  mtext(rev(names_ordered),2,at=1:(dim(plot.dat)[1]),las=2,adj=1,cex=0.85,line = 0.5)
  abline(h=1:(dim(plot.dat)[1])+0.5,lwd=0.3)

  add.legend = F
  if(add.legend){
  image.plot(1:dim(plot.dat)[2],1:dim(plot.dat)[1],round(t(as.matrix(plot.dat)),2),col = rev(mycols(6)),zlim=c(-1,1),breaks = c(-1,-0.65,-0.3,0,0.3,0.65,1),
             horizontal=T, xlab='',legend.only = T,
             axis.args=list(at=tick,
                            labels=c('-1.00','-0.65','-0.30',
                                     '0','0.30','0.65','1.00')),
             legend.cex = 0.8)   
  }
  
  
}


plot.loadings2.old <- function(fa_output, factor.names = NULL){
  # fa_output is the output from fa()
  # factor.names is a character vector that user can supply to rename the factors
  # use 0.3 as cutoff
  loadings <- fa_output$loadings
  dat = as.matrix(data.frame(matrix(as.numeric(loadings), attributes(loadings)$dim, dimnames=attributes(loadings)$dimnames)))
  # name the factors
  if(!is.null(factor.names)){
    colnames(dat) <- factor.names
  }

  require(fields)  
  mycols = colorRampPalette(c(tim.colors(10)[c(1,2)],'white','white',tim.colors(10)[c(9,10)]))
  
  tick<- c(-1,-0.65,-0.3,0,0.3,0.65,1)
  par(mar=c(5,10,8,2.5) + 0.1)
  
  plot.dat = dat[rev(rownames(dat)),]
  
  image(1:10,1:27,round(t(as.matrix(plot.dat)),2),col = rev(mycols(6)),zlim=c(-1,1),
        breaks = c(-1,-0.65,-0.3,0,0.3,0.65,1),
        yaxt='n',xaxt='n',
        xlab='',ylab='',cex.lab=1.2,cex.axis=1.4)
  abline(v=1:10+0.5,lwd=1.5)
  box(lwd=1.5)
  #axis(3, labels = FALSE)
  ## Create some text labels
  labels <- colnames(dat)
  ## Plot x axis labels at default tick marks
  text(x=1:10-0.7,y=rep(28,10), adj = 0, pos = 4,offset = 1, srt = 70, 
       labels = labels, xpd = TRUE,cex=0.8)
  
  #mtext(colnames(dat),3,at=1:10,las=2,adj=0,cex=0.8,line = 0.5,srt = 45)
  mtext(rev(names_ordered),2,at=1:27,las=2,adj=1,cex=0.7,line = 0.5)
  abline(h=1:27+0.5,lwd=0.5)
  
}



plot.legend = F
if(plot.legend){
  #### for legend
  image.plot(1:10,1:27,round(t(as.matrix(dat)),2),col = rev(mycols(8)),zlim=c(-1,1),breaks = seq(-1,1,0.25),
             horizontal=T, xlab='',legend.only = T,
             axis.args=list(at=tick,
                            labels=c('-1.00','-0.75','-0.50',
                                     '-0.25','0','0.25','0.50','0.75','1.00')),
             legend.cex = 0.8)      

}


plot.legend = F
if(plot.legend){
  #### for legend
  require(fields)
  loadings <- factor$loadings
  dat = as.matrix(data.frame(matrix(as.numeric(loadings), attributes(loadings)$dim, dimnames=attributes(loadings)$dimnames)))
  # name the factors
  if(!is.null(factor.names)){
    colnames(dat) <- factor.names
  }
  
  mycols = colorRampPalette(c(tim.colors(10)[c(1,2)],'white','white',tim.colors(10)[c(9,10)]))
  
  tick<- c(-1,-0.65,-0.3,0,0.3,0.65,1)
  par(mar=c(5,10,8,2.8) + 0.1)
  
  plot.dat = dat[rev(rownames(dat)),]
  
  image.plot(1:dim(plot.dat)[2],1:dim(plot.dat)[1],round(t(as.matrix(plot.dat)),2),col = rev(mycols(6)),zlim=c(-1,1),breaks = c(-1,-0.65,-0.3,0,0.3,0.65,1),
             horizontal=T, xlab='',legend.only = T,
             axis.args=list(at=tick,
                            labels=c('-1.00','-0.65','-0.30',
                                     '0','0.30','0.65','1.00')),
             legend.cex = 0.8)      
  
}

reorder <- function(clustering){
  # reorder the cluster labels according to size, for comparing two clustering solutions later
  ord = order(table(clustering))
  return(ord[clustering])
}

clearspace<-function(){
  # remove all the object from global environment
  rm(list = ls(envir = .GlobalEnv),envir = .GlobalEnv)
}

# compute validity metrics for a range of k value
testK <- function(cluster.df, k.range = 2:15, method = 'kmeans'){
  # method can be kmeans or pam
  # returns the inverse dbi and asw index as a function of k
  require(clusterCrit) # for cluster internal and external validity metrics
  require(cluster)
  res = NULL
  for(k in k.range){
    set.seed(10) # so results are replicable each time
    if(method == 'kmeans'){
      reduced <-kmeans(cluster.df, k, iter.max = 1000,nstart=1000,) 
    }
    
    if(method == 'pam'){
      reduced = pam(cluster.df,k,stand = F) # do not standardize, if need to standardize, needs to be applied in the input dataframe before calling the function
      reduced$cluster = reduced$clustering
    }
    
    # excluding singular cluster obs when computing metrics, otherwise asw won't be computed
    freq = table(reduced$cluster)
    sig.cls = as.numeric(names(freq))[freq!=1]
    
    tmp1 = as.data.frame(intCriteria(as.matrix(cluster.df[reduced$cluster %in% sig.cls,]),
                                     reduced$cluster[reduced$cluster %in% sig.cls],"all"))
    tmp1$k = k
    rm(sig.cls)
    res = rbind(res,tmp1)
  }
  
  # get inverse dbi, so that quality maximizes the index
  res$davies_bouldin = 1/res$davies_bouldin
  return(res)
  
}





run.TF = F
if(run.TF){
  nhts.vars = read_excel(path = "../Data/Documentation/NHTS/thsc-nhts17-caltrans-codebook.xlsx", sheet = 'Variables')
  nhts.vars = nhts.vars[,1:3]
  colnames(nhts.vars) = c('name','question.label','question.text')
  nhts.vars$name = tolower(nhts.vars$name)
  
  nhts.vals = read_excel(path = "../Data/Documentation/NHTS/thsc-nhts17-caltrans-codebook.xlsx", sheet = 'Value Lookup')
  colnames(nhts.vals) = c('name','table','value','label')
  nhts.vals[,1:2] = apply(nhts.vals[,1:2], c(1,2), tolower)
  
  lookuptable17 = merge(nhts.vars, nhts.vals, all = T)
  lookuptable17 = lookuptable[,c('name','table','value','label','question.label','question.text')]
  
  save(lookuptable17, file = 'nths.codebook.17.RData')
    
}

nhts.lookup.var<- function(var, table = NULL, year = 2017){
  # function to look up values of a variable in nhts data
  # note that 2009 data has no excel codebook, so there is not searchable dictionary to use
  if(year == 2017){
    load('nths.codebook.17.RData')
    lookuptable = lookuptable17
  }else{
    stop('no year found')
  }
    if(var %in% unique(lookuptable$name)){
      if(is.null(table)){
        return(lookuptable[lookuptable$name == var,])
      }else{
        if(table %in% unique(lookuptable$table)){
          return(lookuptable[lookuptable$name == var & lookuptable$table == table,])
        }else{
          stop('no table found')
        }
      }
    }else{
      stop('no variable found')
    }
}

map_labels <- function(dat,raw_colname = 'o_microtype',microtype = T){
  # function to map the raw cluster labels to relabeled names for micro and geotypes
  # input 
  #   dat is a dataframe that contains the raw cluster label in raw_colname
  #   if microtype = T, then take the raw label and map it to microtype labels, otherwise, map to geotype labels
  # return the updated label vector
  require(tidyr); require(dplyr)
  
  dat$inlabel = dat[,raw_colname]
  if(microtype){
    dat = dat %>%
      mutate(updatelabel = case_when(
        as.numeric(as.character(inlabel)) == 1 ~ as.numeric(4),
        as.numeric(as.character(inlabel)) == 2 ~ as.numeric(2),
        as.numeric(as.character(inlabel)) == 3 ~ as.numeric(6),
        as.numeric(as.character(inlabel)) == 4 ~ as.numeric(3),
        as.numeric(as.character(inlabel)) == 5 ~ as.numeric(5),
        as.numeric(as.character(inlabel)) == 6 ~ as.numeric(1),
        as.numeric(as.character(inlabel)) == 7 ~ as.numeric(0)
      ))
  }else{
    dat = dat %>%
      mutate(updatelabel = case_when(
        as.numeric(as.character(inlabel)) == 1 ~ as.character('G'),
        as.numeric(as.character(inlabel)) == 2 ~ as.character('E'),
        as.numeric(as.character(inlabel)) == 3 ~ as.character('H'),
        as.numeric(as.character(inlabel)) == 4 ~ as.character('D'),
        as.numeric(as.character(inlabel)) == 5 ~ as.character('B'),
        as.numeric(as.character(inlabel)) == 6 ~ as.character('F'),
        as.numeric(as.character(inlabel)) == 7 ~ as.character('C'),
        as.numeric(as.character(inlabel)) == 8 ~ as.character('I'),
        as.numeric(as.character(inlabel)) == 9 ~ as.character('A')
      ))
    
    
  }
  return(dat$updatelabel)
}

 
run7modes <- function(dat){
  for(trip_purp in trip_purps){
    print(trip_purp)
    tmpdat = dat[dat$trip_purp == trip_purp,]
    table(tmpdat$mode, tmpdat$choice)
    tm1 <- mlogit.data(tmpdat, choice = "choice", shape = "long", 
                       chid.var = "trip_indx", alt.var = "mode", drop.index = TRUE)
    
    
    tic()
    ml.tm1 <- mlogit(f1, tm1,weights = wtperfin,
                     nests = list(transit = c("bus", "taxi", "rail_l", "rail_c"), 
                                  auto = c("hv"), micromobility = c("bike", "walk")), un.nest.el = TRUE)
    
    toc() # 
    #show results
    summary(ml.tm1)
    
    # store results
    ml.list[[trip_purp]] <<- ml.tm1
    
    rm(ml.tm1)
    
    
  }
  
}
  
run5modes <- function(dat){
  for(trip_purp in trip_purps){
    print(trip_purp)
    tmpdat = dat[dat$trip_purp == trip_purp,]
    table(tmpdat$mode, tmpdat$choice)
    tm1 <- mlogit.data(tmpdat, choice = "choice", shape = "long", 
                       chid.var = "trip_indx", alt.var = "mode", drop.index = TRUE)
    
    
    tic()
    ml.tm1 <- mlogit(f1.norail, tm1,weights = wtperfin,
                     nests = list(transit = c("bus", "taxi"), 
                                  auto = c("hv"), micromobility = c("bike", "walk")), un.nest.el = TRUE)
    
    toc() # 
    #show results
    summary(ml.tm1)
    
    # store results
    ml.list[[trip_purp]] <<- ml.tm1
    
    rm(ml.tm1)
    
    
  }
  
}


run7modes2 <- function(dat, f){
  # nested logit with formula as input
  for(trip_purp in trip_purps){
    print(trip_purp)
    tmpdat = dat[dat$trip_purp == trip_purp,]
    table(tmpdat$mode, tmpdat$choice)
    tm1 <- mlogit.data(tmpdat, choice = "choice", shape = "long", 
                       chid.var = "trip_indx", alt.var = "mode", drop.index = TRUE)
    
    
    tic()
    tryCatch({
      ml.tm1 <- mlogit(f, tm1,weights = wtperfin,
                       nests = list(transit = c("bus", "taxi", "rail_l", "rail_c"), 
                                    auto = c("hv"), micromobility = c("bike", "walk")), un.nest.el = TRUE)
      
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    
    toc() # 
    #show results
#    summary(ml.tm1)
    
    # store results
    tryCatch({
      
      ml.list[[trip_purp]] <<- ml.tm1
      rm(ml.tm1)
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})

    
  }
  
}


run7modes2.nonesting <- function(dat, f=NULL){
  for(trip_purp in trip_purps){
    print(trip_purp)
    tmpdat = dat[dat$trip_purp == trip_purp,]
    table(tmpdat$mode, tmpdat$choice)
    tm1 <- mlogit.data(tmpdat, choice = "choice", shape = "long", 
                       chid.var = "trip_indx", alt.var = "mode", drop.index = TRUE)
    
    
    tic()
    if(is.null(f)){
      ml.tm1 <- mlogit(choice ~ inv_time + access_time + wait_time |1, tm1,
                       weights = wtperfin, reflevel = "hv")
    }else{
    ml.tm1 <- mlogit(f, tm1,weights = wtperfin, reflevel = "hv")
    }
    
    toc() # 
    #show results
    summary(ml.tm1)
    
    # store results
    ml.list[[trip_purp]] <<- ml.tm1
    
    rm(ml.tm1)
    
    
  }
  
}


### plot the waittime and mode choice contribution, nonlinear relationship. 
plot_waittm_nonlinear <- function(a, aa){
  x = seq(0,20, 0.1)
  y = a*x + aa*x^2
  plot(y~x)
}

#########
testTF = F
if(testTF){
  tdat = dat2 %>% 
    filter(MicrotypeID == 'A_3', f_sys_bins == 'Freeways', fips == 17) %>%
    select(-Tmc, -Miles, - AADT) %>%
    distinct()
  
  density = tdat$census_density
  flow = tdat$census_flow
}


fit.mfd <- function(input.density, input.flow, method = 'three-steps'){
  # input is the data points on flow-density plot
  # fitting method: method = 'three-steps' or 'quadratic'
  # output is a dataframe:
  # res$freeflow.speed = freeflow.speed (the slope)
  # res$freeflow.endpoint = c(density, flow)  # the start point is (0,0)
  # 
  # res$capacity = c( critical density, maximum flow value)
  # congestion phase results are set to NA if there are less than 20 data points observed after that
  # res$gridlock.density = gridlock density (the intercept with x axis)
  
  require(dplyr)
  require(tidyr)
  res = NULL
  
  dat = data.frame(density = input.density, flow = input.flow)
  dat$speed = dat$flow/dat$density

  if(method == 'three-steps'){
    # step 1: fit the free flow curve, a line passing origin:
    # get the 85 percentile of speed
    spd.85 = quantile(dat$speed, probs = 0.85)
    sub = dat %>%
      filter(speed >= spd.85)
    
    tmp = lm(flow~density -1 , data = sub)
    
    res$freeflow.speed = tmp$coefficients
    
    # free flow end point
    ff.maxflow = max(sub$flow)
    ff.maxdensity = max(sub[sub$flow == ff.maxflow,'density'])
    
    res$freeflow.endpoint = c(ff.maxdensity, ff.maxflow)
    
    # step 2: get the capacity: maximum flow value, and critical density
    
    maxflow = max(dat$flow)
    maxdensity = min(dat[dat$flow == maxflow, 'density'])
    
    res$capacity = c(maxdensity, maxflow)
    
    # step 3: get the gridlock density (the intercept with x axis)
    
    rm(tmp)
    rm(sub)
    # get the data points beyond the critical density
    sub = dat %>% filter(density > maxdensity)
    if(dim(sub)[1]<20){
      res$gridlock.density = NA
    }else{
    
      tmp = lm( I(flow-maxflow) ~ I(density-maxdensity) + 0, data = sub)
      
      res$gridlock.density = maxdensity - maxflow/tmp$coefficients
      
    #  plot(dat$density, dat$flow, type = 'b', xlim = c(0, 1000))
    #  lines(sub$density, tmp$coefficients*(sub$density - maxdensity)+maxflow, col='blue')
    
      testing = F
      if(testing){
        ggplot(data = dat, aes(x = density, y = flow)) +
        geom_point() +
        geom_segment(aes(x= c(0), y = c(0), 
                         xend = c(res$freeflow.endpoint[1]),yend =c(res$freeflow.endpoint[2])), 
                     color = 'green' , size = rel(1.4)) +
        geom_segment(aes(x= c(res$freeflow.endpoint[1]), y = c(res$freeflow.endpoint[2]), 
                       xend = c( res$capacity[1]),yend =c( res$capacity[2])), 
                   color = 'green' , size = rel(1.4)) +
        geom_segment(aes(x = c( res$capacity[1]),y =c( res$capacity[2]), 
                         xend = c(res$gridlock.density),yend =c(0 )), 
                     color = 'pink' , size = rel(1.4)) +
          coord_cartesian(xlim=c(0, max(dat$density)))
      }
    }
    return.res = data.frame(ff.spd = res$freeflow.speed, 
                            ff.maxdensity, ff.maxflow, 
                            capacity.flow = maxflow,
                            capacity.density = maxdensity,
                            gridlock.density = res$gridlock.density)
  }
  
  if(method == 'quadratic'){
    # now fit a flow = phi*density(density - gridlock.density)
    
    return.res = lm(data = dat, flow ~ I(density^2) + density -1)
    
  }
  
  if(method == 'polynomial'){
    # now fit a flow = phi*density(density - gridlock.density)
    
    return.res = lm(data = dat, flow ~ I(density^3) + I(density^2) + density -1)
    
  }
  
  
  return(return.res)
  
}

fit.mfd.rural <- function(x, y){
  # input is the data points on flow-density plot
  # x = density, y = flow
  # we find out freeflow speed and capacity flow value and return these two parameters. 
  # output is a dataframe:
  # res$freeflow.speed = freeflow.speed (the slope)
  # res$capacity.flow = capacity flow
  # res$capacity.density = density at the capacity flow.
 
  require(dplyr)
  require(tidyr)

  dat = data.frame(density = x, flow = y)
  dat$speed = dat$flow/dat$density
  
    # step 1: fit the free flow curve, a line passing origin:
    # get the 75 percentile of speed # used to be 85%
  # change to 0.6
    spd = quantile(dat$speed, probs = 0.5)
    sub = dat %>%
      filter(speed >= spd)
    
    tmp = lm(flow~density -1 , data = sub)
    
    ff.spd = tmp$coefficients
    

    # step 2: get the capacity: maximum flow value, and critical density
    
    maxflow = quantile(dat$flow,0.95) # change from 0.99 to 0.95
    maxdensity = maxflow/ff.spd
    

   return.res = data.frame(ff.spd = ff.spd, 
                            capacity.flow = maxflow,
                            capacity.density = maxdensity)
  
  
  return(return.res)
  
}

fit.mfd.rural2 <- function(x, y){
  # input is the data points on flow-density plot
  # x = density, y = flow
  # first find out the capacity flow, then fit the free flow speed 
  # output is a dataframe:
  # res$freeflow.speed = freeflow.speed (the slope)
  # res$capacity.flow = capacity flow
  # res$capacity.density = density at the capacity flow.
  
  require(dplyr)
  require(tidyr)
  
  dat = data.frame(density = x, flow = y)
  dat$speed = dat$flow/dat$density

  
  # step 1: get the capacity: maximum flow value, and critical density
  
  maxflow = quantile(dat$flow,0.95) # change from 0.99 to 0.95
  maxdensity = min(dat$density[dat$flow>=maxflow])
  
  # step 2: fit the free flow curve to the points left of the (maxdensity,maxflow), a line passing origin:
  sub = dat %>%
    filter(density <= maxdensity)
  
  tmp = lm(flow~density -1 , data = sub)
  
  ff.spd = tmp$coefficients
  
  

  return.res = data.frame(ff.spd = ff.spd, 
                          capacity.flow = maxflow,
                          capacity.density = maxdensity)
  
  
  return(return.res)
  
}

fit.mfd.unified <- function(input.density, input.flow, mfd.shape = 'rural'){
  # the input unit should be density [=] veh/m-lane; flow[=] veh/sec-lane
  # input is the data points on flow-density plot
  # fitting method: mfd.shape = 
  #            'urban': for urban mfd, polynomial fit of order 3 if sufficient data; if not, use triangle fit
  #            'rural': for rural mfd, use step linear
  #            '
  
  require(dplyr)
  require(tidyr)
  res = NULL
  
  dat = data.frame(density = input.density, flow = input.flow)
  dat$speed = dat$flow/dat$density
  
  # step 0: data preprocessing to filter the outter bound
  # done outside of the function before inputs
  
  if(length(input.density)<100){return(NA)}
  # step 1: get the capacity: maximum flow value, and critical density
  
  maxflow = quantile(dat$flow,0.99) # change from 0.99 to 0.95 # change it back to 0.99 because in some cases it would bias the speed high
  maxdensity = min(dat$density[dat$flow>=maxflow])
  
  # step 2: fit the free flow curve to the points left of the (maxdensity,maxflow), a line passing origin:
  sub = dat %>%
    filter(density <= maxdensity)
  
  tmp = lm(flow~density -1 , data = sub)
  
  ff.spd = tmp$coefficients
  
  info.tmp = summary(tmp)
  R2.fit = info.tmp$r.squared 
  
  
  if(mfd.shape == 'rural'){
    # use step linear function
    return.res = data.frame(ff.spd.steplinear = ff.spd, 
                            capacity.flow.steplinear = maxflow,
                            capacity.density.steplinear = maxdensity,
                            R2.fit.steplinear = R2.fit)
    
    
    return(return.res)
    
  }
  
  if(mfd.shape == 'urban'){
# for now fit a parabola
    
    res = lm(data = dat, flow ~ I(density^2) + density -1)
    info.tmp = summary(res)
    R2.fit = info.tmp$r.squared 
    
    # if in the form of ax(x-b), then 
    a = res$coefficients[1]
    b = - res$coefficients[2]/res$coefficients[1]
    
    return(data.frame(coef.x2.parabolic = res$coefficients[1], 
                      coef.x1.parabolic = res$coefficients[2],
                      a.parabolic = a,
                      b.parabolic = b,
                      max.spd.parabolic = res$coefficients[2],
                      jamd.fit.parabolic = - res$coefficients[2]/res$coefficients[1],
                      capcity.fit.parabolic = -a*b^2/4,
                      critical.density.fit.parabolic = -res$coefficients[2]/res$coefficients[1]/2,
                      R2.fit.parabolic = R2.fit))
    
    
  }
  
  if(method == 'polynomial'){
    # now fit a flow = phi*density(density - gridlock.density)
    
    return.res = lm(data = dat, flow ~ I(density^3) + I(density^2) + density -1)
    
  }
  
  
}
  



print.paste <- function(x,...){
  # replacing print(paste(x...))
  print(paste(x,...))
}

# function to further clean the key
keyCleaning <- function(){
  # key is a global variable already read in and joined to transgeo clusters
  dim(key)
  key <<- distinct(key)
  
  print(paste('Raw key has',dim(key)[1],'rows.'))
  dim(key) # [1] 55481    21
  
  key <<- key %>% mutate(dir = case_when(
    dir == 1 ~ 'N',
    dir == 2 ~ 'NE',
    dir == 3 ~ 'E',
    dir == 4 ~ 'SE',
    dir == 5 ~ 'S',
    dir == 6 ~ 'SW',
    dir == 7 ~ 'W',
    dir == 8 ~ 'NW',
    dir == 9 ~ 'N-S or NE-SW',
    dir == 0 ~ 'E-W or SE-NW'))
  
  # combines functional types: road1 = 1&2, road2 = 3&4&5, road3 = 6&7
  key <<- key %>% mutate(f_sys_bins = case_when(
    f_systm %in% c(1,2) ~ 'road1',
    f_systm %in% c(3,4,5) ~ 'road2',
    f_systm %in% c(6,7) ~ 'road3'))
  
  key <<- key %>% mutate(F_sys_bins = case_when(
    F_System %in% c(1,2) ~ 'road1',
    F_System %in% c(3,4,5) ~ 'road2',
    F_System %in% c(6,7) ~ 'road3'))
  
  # how many v stations in the raw key now
  print(paste('Without considering direction, there are', 
        dim(key %>% select(fips, station)%>%distinct())[1], 'vol stns in raw key.'))
  
  
  # do we still have stations fall in more than 1 ftype? only 1.7%, much better.
  tt =key%>% select(fips, station,f_systm) %>% distinct() %>% group_by(fips,station)%>%mutate(nftype = n_distinct(f_systm))
  
  print(paste('There are', sum(tt$nftype>1)/length(tt$nftype), 'vol stns have more than 1 ftype associated with it.'))
  
  # only keep the rows when the direction (dir from volume) matches between volume and speed.
  # and functional types matches
  
  key.clean = key %>%
    filter(dir == Direction) 
  
  print('# further match direction -------------')
  print(paste('There are', dim(key.clean)[1], 'rows left.'))
  print(paste('There are', dim(key.clean %>% select(fips, station)%>%distinct())[1], 'vol stns left'))
  
  key.clean = key.clean %>% filter(F_sys_bins == f_sys_bins) # 19046     22
  print('# further match coarsened functional types -------------')
  print(paste('There are', dim(key.clean)[1], 'rows left.'))
  print(paste('There are', dim(key.clean %>% select(fips, station)%>%distinct())[1], 'vol stns left'))
  
  # remove duplicates
  key.clean = key.clean %>%
    select(-lane) %>%
    distinct()
  print('# remove vol lane index as they associated with their respective vol stations')
  print(paste('There are', dim(key.clean)[1], 'rows left.'))
  print(paste('There are', dim(key.clean %>% select(fips, station)%>%distinct())[1], 'vol stns left'))
  
  dim(key.clean) #9098   22
  
  print('# join to xwalk --------------------')
  key <<- key.clean %>%
    left_join(xwalk)
  print(paste('There are', dim(key)[1], 'rows left.'))
  print(paste('There are', dim(key %>% select(fips, station)%>%distinct())[1], 'vol stns left'))
  
}

merge_cattlab <- function(fname){
  # merge microtype and blobid to cattlab geoids 
  cattlab1 = fread(file.path(datadir, 'CATTLab/MD_PA',paste0(fname,'.csv')))
  
  
  ## merge cattlab
  cattlab1 = cattlab1 %>%
    left_join(transgeo)
  
  
  # select microtype 1, 2
  # cattlab = cattlab %>%
  #   filter(transgeo_microtype %in% c(1,2))
  
  cattlab1 = cattlab1 %>%
    select(-spatial_id )
  
  cattlab1 = cattlab1 %>%
    rename(microtype = transgeo_microtype,
           blobID = transgeo_blobID,
           geotype = transgeo_geotype)
  
  cattlab1 = cattlab1 %>%
    select(-cbsa, -geotype)
  
  return(cattlab1)
  
} 


# function to process data from blob level data from v1 and v2
mfd.proc.v12 <- function( filename, spatial_unit = 'blob'){
  # filename: contains the path and filename
  # 
  dat = fread(filename)
  cnt_tmcs = mean(dat$cnt_tmcs,na.rm = T)
  master.dat = dat
  
  dat1 = master.dat %>% select(-contains("V2"))%>% mutate(version = 'V1') %>%
    rename(agg_flow_plph = V1_aggflow_plph,
           agg_density_plpm = V1_aggdensity_plpm,
           agg_speed_mph = V1_aggspeed_mph)
  dat2 = master.dat %>% select(-contains("V1"))%>% mutate(version = 'V2') %>%
    rename(agg_flow_plph = V2_aggflow_plph,
           agg_density_plpm = V2_aggdensity_plpm,
           agg_speed_mph = V2_aggspeed_mph)
  
  master.dat = bind_rows(dat1,dat2)
  
  if(spatial_unit == 'blob'){
    
    
    cdat = master.dat %>%
      group_by(blobid, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(blobid,microtype, dens.bins) %>%
      mutate(c.flow.l = quantile(agg_flow_plph,0.6),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
  }
  
  if(spatial_unit == 'tract'){
    
    
    cdat = master.dat %>%
      group_by(tract_geoid, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(tract_geoid,microtype, dens.bins) %>%
      mutate(c.flow.l = quantile(agg_flow_plph,0.6),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(tract_geoid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(tract_geoid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
  }
  
  
  res = mfd.res_rural %>%
    left_join(mfd.res_urban)
  res$cnt_tmcs = cnt_tmcs
  return(res)
  
}

# urban mfd function
urban.mfd.curve = Vectorize(function(a, b, d, x){
  #if(x*27.7 < 0.61){return(x*28)}else{return(0.61)}
  if(x>=0 & x<=d-sqrt(d^2-b*d)){return(a*x*(x-b))}else if(x>d-sqrt(d^2-b*d) & x<=d){
    return(a*(2*d - b - 2*sqrt(d^2-b*d))*(x-d))
  }else{return(NA)}
})
# test

rural.mfd.curve = Vectorize(function(a, b, x){
  #curve-linear: constant production
  if(x<b/2){return(a*x*(x-b))}else if(x>=b/2){
    return(-a*b^2/4)
  }else{return(NA)}
})

rural.mfd.linear = Vectorize(function(spd, capacity, x){
  #curve-linear: constant production
  if(x*spd < capacity){return(x*spd)}else{return(capacity)}
})

# microtype 3 parameters
df3.curve = data.frame(x = seq(0,0.15,0.001), y = rural.mfd.curve(-336,0.0658,x=seq(0,0.15,0.001)))
df3.linear = data.frame(x = seq(0,0.15,0.001), y = rural.mfd.linear(spd = 14.4, capacity = 0.354, x = seq(0,0.15,0.001)))
df3.curve$form = 'rural_curve'
df3.linear$form = 'rural_linear'
df3.urban = data.frame(x = seq(0,0.15,0.001), y = urban.mfd.curve(-336,0.0658, 0.2, x=seq(0,0.15,0.001)))
df3.urban$form = 'urban'
df = bind_rows(df3.curve,df3.linear) %>% bind_rows(df3.urban)

ggplot(df, aes(x = x, y = y, color = form))+
  geom_point()

### function to estimate mfd parameters for V4 (80 percentile) and V5 (90 percentile)

mfd.proc.v1245 <- function( filename, spatial_unit = 'blob'){
  # filename: contains the path and filename
  # 
  dat = fread(filename)
  cnt_tmcs = mean(dat$cnt_tmcs,na.rm = T)
  master.dat = dat
  
  dat1 = master.dat %>% select(blobid:datetime,contains("V1") )%>% mutate(version = 'V1') %>%
    rename(agg_flow_plph = V1_aggflow_plph,
           agg_density_plpm = V1_aggdensity_plpm)
  dat2 = master.dat %>% select(blobid:datetime,contains("V2") )%>% mutate(version = 'V2') %>%
    rename(agg_flow_plph = V2_aggflow_plph,
           agg_density_plpm = V2_aggdensity_plpm)
  dat4 = master.dat %>% select(blobid:datetime,contains("V4") )%>% mutate(version = 'V4') %>%
    rename(agg_flow_plph = V4_aggflow_plph,
           agg_density_plpm = V4_aggdensity_plpm)
  dat5 = master.dat %>% select(blobid:datetime,contains("V5") )%>% mutate(version = 'V5') %>%
    rename(agg_flow_plph = V5_aggflow_plph,
           agg_density_plpm = V5_aggdensity_plpm)
  
  master.dat = bind_rows(dat1,dat2, dat4, dat5)
  
  if(spatial_unit == 'blob'){
    
    
    cdat = master.dat %>%
      group_by(blobid, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(blobid,microtype, dens.bins, version) %>%
      mutate(c.flow.l = quantile(agg_flow_plph,0.6),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
  }
  
  if(spatial_unit == 'tract'){
    
    
    cdat = master.dat %>%
      group_by(tract_geoid, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(tract_geoid,microtype, dens.bins) %>%
      mutate(c.flow.l = quantile(agg_flow_plph,0.6),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(tract_geoid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(tract_geoid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
  }
  
  
  res = mfd.res_rural %>%
    left_join(mfd.res_urban)
  res$cnt_tmcs = cnt_tmcs
  return(res)
  
}

mfd.proc.v15 <- function( filename, spatial_unit = 'blob', upperbound_pct = 0.4, zipTF = F, zipfolder, fnm.inzip){
  # filename: contains the path and filename in unzipped folder
  # 
  if(zipTF){
    # if read from a zip file
    library(pacman)
    p_load(R.utils) # for reading gz files
    p_load(readr) # for read files from a zip foler
    
    dat <- read_csv(unzip(zipfolder, fnm.inzip))
    
    
  }else{dat = fread(filename)}
  
  # filter out rows with invalid spd due to unavailable flow density
  dat = dat %>% filter(!is.na(V1_aggspeed_mph))
  
  
  if(nrow(dat)>0){
  cnt_tmcs = mean(dat$cnt_tmcs,na.rm = T)
  master.dat = dat
  
  
  if(spatial_unit == 'blob'){
    
    dat1 = master.dat %>% select(blob_id:datetime,contains("V1") )%>% mutate(version = 'V1') %>%
      rename(agg_flow_plph = V1_aggflow_plph,
             agg_density_plpm = V1_aggdensity_plpm)
    # dat2 = master.dat %>% select(blobid:datetime,contains("V2") )%>% mutate(version = 'V2') %>%
    #   rename(agg_flow_plph = V2_aggflow_plph,
    #          agg_density_plpm = V2_aggdensity_plpm)
    # dat4 = master.dat %>% select(blobid:datetime,contains("V4") )%>% mutate(version = 'V4') %>%
    #   rename(agg_flow_plph = V4_aggflow_plph,
    #          agg_density_plpm = V4_aggdensity_plpm)
    dat5 = master.dat %>% select(blob_id:datetime,contains("V5") )%>% mutate(version = 'V5') %>%
      rename(agg_flow_plph = V5_aggflow_plph,
             agg_density_plpm = V5_aggdensity_plpm)
    
    master.dat = bind_rows(dat1, dat5)
    
    cdat = master.dat %>%
      group_by(blob_id, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(blob_id,microtype, dens.bins, version) %>%
      mutate(c.flow.l = quantile(agg_flow_plph,(1-upperbound_pct)),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(blob_id, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(blob_id, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
  }
  
  if(spatial_unit == 'tract'){
    
    dat1 = master.dat %>% select(tract_geoid:datetime,contains("V1") )%>% mutate(version = 'V1') %>%
      rename(agg_flow_plph = V1_aggflow_plph,
             agg_density_plpm = V1_aggdensity_plpm)
    # dat2 = master.dat %>% select(blobid:datetime,contains("V2") )%>% mutate(version = 'V2') %>%
    #   rename(agg_flow_plph = V2_aggflow_plph,
    #          agg_density_plpm = V2_aggdensity_plpm)
    # dat4 = master.dat %>% select(blobid:datetime,contains("V4") )%>% mutate(version = 'V4') %>%
    #   rename(agg_flow_plph = V4_aggflow_plph,
    #          agg_density_plpm = V4_aggdensity_plpm)
    dat5 = master.dat %>% select(tract_geoid:datetime,contains("V5") )%>% mutate(version = 'V5') %>%
      rename(agg_flow_plph = V5_aggflow_plph,
             agg_density_plpm = V5_aggdensity_plpm)
    
    master.dat = bind_rows(dat1, dat5)
    
    cdat = master.dat %>%
      group_by(tract_geoid, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(tract_geoid,microtype, dens.bins, version) %>%
      mutate(c.flow.l = quantile(agg_flow_plph, (1-upperbound_pct)),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(tract_geoid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(tract_geoid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
  }
  
  
  res = mfd.res_rural %>%
    left_join(mfd.res_urban)
  res$cnt_tmcs = cnt_tmcs
  return(res)
  } # if the file contains data
  
}


mfd.proc.extractv1 <- function( filename, spatial_unit = 'blob', upperbound_pct = 0.4, zipTF = F, zipfolder, fnm.inzip){
  # extract v1 data in upper 20percentile of each density bin

  if(zipTF){
    # if read from a zip file
    library(pacman)
    p_load(R.utils) # for reading gz files
    p_load(readr) # for read files from a zip foler
    
    dat <- read_csv(unzip(zipfolder, fnm.inzip))
    
    
  }else{dat = fread(filename)}
  
  # filter out rows with invalid spd due to unavailable flow density
  dat = dat %>% filter(!is.na(V1_aggspeed_mph) | V1_aggdensity_plpm < 0.00001)
  
  
  if(nrow(dat)>0){
    cnt_tmcs = mean(dat$cnt_tmcs,na.rm = T)
    master.dat = dat
    
    
    if(spatial_unit == 'blob'){
      
      dat1 = master.dat %>% select(blob_id:datetime,contains("V1") )%>% mutate(version = 'V1') %>%
        rename(agg_flow_plph = V1_aggflow_plph,
               agg_density_plpm = V1_aggdensity_plpm,
               aggspeed_mph = V1_aggspeed_mph )
      
      

      master.dat = bind_rows(dat1)
      # error handling of bad data
      if(median(master.dat$agg_density_plpm,na.rm = T)<=0.001 | median(master.dat$agg_flow_plph,na.rm = T) <0.001){
        return(NULL)
      }
      
      cdat = master.dat %>%
        group_by(blob_id, microtype, version) %>%
        mutate(nsample = n())%>%
        filter(nsample>1000)%>%
        select(-nsample)%>%
        mutate(maxdens = max(agg_density_plpm,na.rm = T),
               mindens = min(agg_density_plpm, na.rm = T))%>%
        mutate( dens.bins = cut(agg_density_plpm, 
                                breaks = seq(unique(mindens),unique(maxdens),
                                             by = (unique(maxdens)-unique(mindens))/60))) %>%
        ungroup()%>%
        group_by(blob_id,microtype, dens.bins, version) %>%
        mutate(c.flow.l = quantile(agg_flow_plph,(1-upperbound_pct)),
               c.flow.h = quantile(agg_flow_plph,0.99),
               nsample = n()) %>%
        ungroup() %>%
        filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
        dplyr::select(-c.flow.l, -c.flow.h, -nsample, -maxdens, -mindens,-version, -dens.bins)
   #     mutate(agg_flow_lane_sec = agg_flow_plph /3600,
  #             agg_density_lane_meter = agg_density_plpm/1609.34) # only extract data, no need to convert unit
   }
    
    if(spatial_unit == 'tract'){
      
      dat1 = master.dat %>% select(tract_geoid:datetime,contains("V1") )%>% 
        mutate(version = 'V1') %>%
        rename(agg_flow_plph = V1_aggflow_plph,
               agg_density_plpm = V1_aggdensity_plpm)
 
      master.dat = bind_rows(dat1)
      
      # error handling of bad data
      if(median(master.dat$agg_density_plpm,na.rm = T)<=0.001 | median(master.dat$agg_flow_plph,na.rm = T) <0.001){
        return(NULL)
      }
      
      cdat = master.dat %>%
        group_by(tract_geoid, microtype, version) %>%
        mutate(nsample = n())
      if(cdat$nsample[1]<1000)return(NULL)
      cdat = cdat %>%
        filter(nsample>1000)%>%
        select(-nsample)%>%
        mutate(maxdens = max(agg_density_plpm,na.rm = T),
               mindens = min(agg_density_plpm, na.rm = T))%>%
        mutate( dens.bins = cut(agg_density_plpm, 
                                breaks = seq(unique(mindens),unique(maxdens),
                                             by = (unique(maxdens)-unique(mindens))/60))) %>%
        ungroup()%>%
        group_by(tract_geoid,microtype, dens.bins, version) %>%
        mutate(c.flow.l = quantile(agg_flow_plph, (1-upperbound_pct)),
               c.flow.h = quantile(agg_flow_plph,0.99),
               nsample = n()) %>%
        ungroup() %>%
        filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
        dplyr::select(-c.flow.l, -c.flow.h, -nsample, -maxdens, -mindens,-version, -dens.bins)
      
     }
    
    cdat$cnt_tmcs = cnt_tmcs
    return(cdat)
  } # if the file contains data
  
}



# make plot of raw and upper bound and overlay the fit
mfd.plot.v1 <- function( filename, spatial_unit = 'blob', upperbound_pct = 0.4){
  # filename: contains the path and filename
  # 
  dat = fread(filename)
  sblobid = unique(dat$blobid)
  xwalk.tmp = xwalk.blob %>% filter(blobid == sblobid)
  sgeotype = xwalk.tmp$geotype[1]
  smicrotype = xwalk.tmp$microtype[1]
  
  cnt_tmcs = mean(dat$cnt_tmcs,na.rm = T)
  master.dat = dat
  
  dat1 = master.dat %>% select(blobid:datetime,contains("V1") )%>% mutate(version = 'V1') %>%
    rename(agg_flow_plph = V1_aggflow_plph,
           agg_density_plpm = V1_aggdensity_plpm)
  dat1$version = 'V1'
  master.dat = dat1
  
  dat1 = dat1 %>%
    mutate(y = agg_flow_plph /3600,
           x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
  

  if(spatial_unit == 'blob'){
    
    
    cdat = master.dat %>%
      group_by(blobid, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(blobid,microtype, dens.bins, version) %>%
      mutate(c.flow.l = quantile(agg_flow_plph,(1-upperbound_pct)),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
    
    # make plots of raw and upperbound scatter, overlay with urban curves
    
    myfun = function(x){
      #if(x*27.7 < 0.61){return(x*28)}else{return(0.61)}
      return(mfd.res_urban$a.parabolic*x^2 - mfd.res_urban$a.parabolic*mfd.res_urban$b.parabolic*x)
    }
    
    ggplot(cdat, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.16)) +
      scale_y_continuous(limits = c(0, 0.55)) +
      theme_bw()+
      
      ggtitle(paste('upperbound',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('upperbd',sgeotype,smicrotype,'blob',sblobid,'quadform.pdf',sep = '.')),height = 3, width = 3.5)
 
    ggplot(dat1, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.16)) +
      scale_y_continuous(limits = c(0, 0.55)) +
      theme_bw()+
      
      ggtitle(paste('rawdat',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('rawdat',sgeotype,smicrotype,'blob',sblobid,'quadform.pdf',sep = '.')),height = 3, width = 3.5)
    
    
    # quatratic + backward wave
    myfun = Vectorize(function(x){
      #if(x*27.7 < 0.61){return(x*28)}else{return(0.61)}
      if(x<=0.15-sqrt(0.15^2-mfd.res_urban$b.parabolic*0.15)){
      return(mfd.res_urban$a.parabolic*x^2 - mfd.res_urban$a.parabolic*mfd.res_urban$b.parabolic*x)}
      else{
        return(mfd.res_urban$a.parabolic * (2*0.15 - mfd.res_urban$b.parabolic - 2*sqrt(0.15^2-mfd.res_urban$b.parabolic*0.15))*(x-0.15))
      }
    })
    
    ggplot(cdat, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.16)) +
      scale_y_continuous(limits = c(0, 0.55)) +
      theme_bw()+
      
      ggtitle(paste('upperbound',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('upperbd',sgeotype,smicrotype,'blob',sblobid,'quad_urban.pdf',sep = '.')),height = 3, width = 3.5)

    
    # quatratic + fixed capacity
    myfun = Vectorize(function(x){
      #if(x*27.7 < 0.61){return(x*28)}else{return(0.61)}
      if(x<= mfd.res_urban$b.parabolic/2){
        return(mfd.res_urban$a.parabolic*x^2 - mfd.res_urban$a.parabolic*mfd.res_urban$b.parabolic*x)}
      else{
        return(-mfd.res_urban$a.parabolic *mfd.res_urban$b.parabolic^2/4)
      }
    })
    
    ggplot(cdat, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.16)) +
      scale_y_continuous(limits = c(0, 0.55)) +
      theme_bw()+
      
      ggtitle(paste('upperbound',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('upperbd',sgeotype,smicrotype,'blob',sblobid,'quad_rural.pdf',sep = '.')),height = 3, width = 3.5)
    
      
      
    
  }
  
}

mfd.plot.v1.ppt <- function( filename, spatial_unit = 'blob', upperbound_pct = 0.4){
  # filename: contains the path and filename
  # 
  dat = fread(filename)
  sblobid = unique(dat$blobid)
  xwalk.tmp = xwalk.blob %>% filter(blobid == sblobid)
  sgeotype = xwalk.tmp$geotype[1]
  smicrotype = xwalk.tmp$microtype[1]
  
  cnt_tmcs = mean(dat$cnt_tmcs,na.rm = T)
  master.dat = dat
  
  dat1 = master.dat %>% select(blobid:datetime,contains("V1") )%>% mutate(version = 'V1') %>%
    rename(agg_flow_plph = V1_aggflow_plph,
           agg_density_plpm = V1_aggdensity_plpm)
  dat1$version = 'V1'
  master.dat = dat1
  
  dat1 = dat1 %>%
    mutate(y = agg_flow_plph /3600,
           x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
  
  
  if(spatial_unit == 'blob'){
    
    
    cdat = master.dat %>%
      group_by(blobid, microtype, version) %>%
      mutate(nsample = n())%>%
      filter(nsample>1000)%>%
      select(-nsample)%>%
      mutate(maxdens = max(agg_density_plpm,na.rm = T),
             mindens = min(agg_density_plpm, na.rm = T))%>%
      mutate( dens.bins = cut(agg_density_plpm, 
                              breaks = seq(unique(mindens),unique(maxdens),
                                           by = (unique(maxdens)-unique(mindens))/60))) %>%
      ungroup()%>%
      group_by(blobid,microtype, dens.bins, version) %>%
      mutate(c.flow.l = quantile(agg_flow_plph,(1-upperbound_pct)),
             c.flow.h = quantile(agg_flow_plph,0.99),
             nsample = n()) %>%
      ungroup() %>%
      filter((agg_flow_plph >= c.flow.l & agg_flow_plph <= c.flow.h ) | nsample <=3 ) %>%
      mutate(y = agg_flow_plph /3600,
             x = agg_density_plpm/1609.34) # convert to veh/m-lane, and veh/s-lane
    
    mfd.res_rural = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'rural'))
    
    mfd.res_urban = cdat %>%
      group_by(blobid, microtype,version)%>%
      do(fit.mfd.unified(input.density = .$x,input.flow = .$y, mfd.shape = 'urban'))
    
    # make plots of raw and upperbound scatter, overlay with urban curves
    
    myfun = function(x){
      #if(x*27.7 < 0.61){return(x*28)}else{return(0.61)}
      return(mfd.res_urban$a.parabolic*x^2 - mfd.res_urban$a.parabolic*mfd.res_urban$b.parabolic*x)
    }
    
    ggplot(cdat, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.2)) +
      scale_y_continuous(limits = c(0, 0.45)) +
      theme_bw()+
      
      ggtitle(paste('upperbound',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('ppt.upperbd',sgeotype,smicrotype,'blob',sblobid,'quadform.pdf',sep = '.')),height = 2.5, width = 4.5)
    
    ggplot(dat1, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.2)) +
      scale_y_continuous(limits = c(0, 0.45)) +
      theme_bw()+
      
      ggtitle(paste('rawdat',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('ppt.rawdat',sgeotype,smicrotype,'blob',sblobid,'quadform.pdf',sep = '.')),height = 2.5, width = 4.5)
    
    
    # quatratic + backward wave
    myfun = Vectorize(function(x){
      #if(x*27.7 < 0.61){return(x*28)}else{return(0.61)}
      if(x<=0.15-sqrt(0.15^2-mfd.res_urban$b.parabolic*0.15)){
        return(mfd.res_urban$a.parabolic*x^2 - mfd.res_urban$a.parabolic*mfd.res_urban$b.parabolic*x)}
      else{
        return(mfd.res_urban$a.parabolic * (2*0.15 - mfd.res_urban$b.parabolic - 2*sqrt(0.15^2-mfd.res_urban$b.parabolic*0.15))*(x-0.15))
      }
    })
    
    ggplot(cdat, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.2)) +
      scale_y_continuous(limits = c(0, 0.45)) +
      theme_bw()+
      
      ggtitle(paste('upperbound',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('ppt.upperbd',sgeotype,smicrotype,'blob',sblobid,'quad_urban.pdf',sep = '.')),height = 2.5, width = 4.5)
    
    
    # quatratic + fixed capacity
    myfun = Vectorize(function(x){
      #if(x*27.7 < 0.61){return(x*28)}else{return(0.61)}
      if(x<= mfd.res_urban$b.parabolic/2){
        return(mfd.res_urban$a.parabolic*x^2 - mfd.res_urban$a.parabolic*mfd.res_urban$b.parabolic*x)}
      else{
        return(-mfd.res_urban$a.parabolic *mfd.res_urban$b.parabolic^2/4)
      }
    })
    
    ggplot(cdat, aes(x = x, y = y)) +
      geom_point() +
      stat_function(fun = myfun, colour = 'red', size = rel(2)) +
      xlab('density (veh/m-lane)')+
      ylab('flow (veh/s-lane)') +
      scale_x_continuous(limits = c(0, 0.2)) +
      scale_y_continuous(limits = c(0, 0.45)) +
      theme_bw()+
      
      ggtitle(paste('upperbound',sgeotype,smicrotype,'blob',sblobid,'R2=',round(mfd.res_urban$R2.fit.parabolic,4) ))
    ggsave(file = file.path(figuredir,paste('ppt.upperbd',sgeotype,smicrotype,'blob',sblobid,'quad_rural.pdf',sep = '.')),height = 2.5, width = 4.5)
    
    
    
    
  }
  
}

emfd.fit <- function(spd, emis){
  require(tidyverse)
  df = data.frame(x = spd, y = emis)
  df = df%>%
    filter(x>0)%>%
    mutate(rx = 1/x)
  
  # fit = nls(data = df , y ~ a/x + b, 
  #           start = list(a = 1, b = 0.2))
  # return(data.frame(a = coef(fit)[1], b = coef(fit)[2]))
  
  # use linear fit now to get R2 reported:
  
  fit = lm(data = df, y ~ rx)
  return(data.frame(a = as.numeric(coef(fit)[2]), b = as.numeric(coef(fit)[1]), adj.r2 = as.numeric(summary(fit)$adj.r.squared)))
  
  
  
}

# emfd.proc <- function(tractid, dat){
#   # tractid is the tract geoid
#   # dat is the processed emfd data
#   tmp <- dat %>% filter(tract_geoid == tractid)
#   res = tmp %>%
#     mutate(spd_km_h=V1_aggspeed_mph*1.609, emis_g_veh_meter = emrate/1609)%>%
#     group_by(modelyear, tract_geoid, microtype, geotype)%>%
#     do(emfd.fit(spd = .$spd_km_h,emis = .$emis_g_veh_meter))
#   
#   return(res)
# }


# function to do ANOVA analysis

gems.anova <- function(dat, dependent.var, smape.TF = F){
   res = data.frame(metrics = c('adj.r.squared','p.value'), microtype = NA, geotype = NA, geomicrotype = NA, geomicrotype.state = NA, geomicrotype.state.smape = NA)
  # 1. microtype alone
  f = as.formula(paste(dependent.var, 'as.factor(microtype)',sep = '~'))
  
  lm1 = lm(f, data = dat,weights = cnt_tmcs)
  res$microtype[res$metrics=='adj.r.squared'] = summary(lm1)$adj.r.squared
  
  # 2. geotype alone
  f = as.formula(paste(dependent.var, 'as.factor(geotype)',sep = '~'))
  
  lm2 = lm(f, data = dat,weights = cnt_tmcs)
  res$geotype[res$metrics=='adj.r.squared'] = summary(lm2)$adj.r.squared
  
  # 3. microtype + geotype
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype)',sep = '~'))
  
  lm3 = lm(f, data = dat,weights = cnt_tmcs)
  res$geomicrotype[res$metrics=='adj.r.squared'] = summary(lm3)$adj.r.squared
  
  # anova 3 to 1
  res$geomicrotype[res$metrics=='p.value']= anova(lm1,lm3)[['Pr(>F)']][2]
  

  # 4. microtype + geotype + state
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype) + as.factor(state)',sep = '~'))
  
  lm4 = lm(f, data = dat,weights = cnt_tmcs)
  res$geomicrotype.state[res$metrics=='adj.r.squared'] = summary(lm4)$adj.r.squared
  
  # anova 4 to 3
  res$geomicrotype.state[res$metrics=='p.value']= anova(lm3,lm4)[['Pr(>F)']][2]
  
  # 5. microtype + geotype + state + smape by frcs
  if(smape.TF){
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype) + as.factor(state) + frc1 + frc2 + frc3 + frc4p',sep = '~'))
  
  lm5 = lm(f, data = dat,weights = cnt_tmcs)
  res$geomicrotype.state.smape[res$metrics=='adj.r.squared'] = summary(lm5)$adj.r.squared
  
  # anova 5to 4
  res$geomicrotype.state.smape[res$metrics=='p.value']= anova(lm4,lm5)[['Pr(>F)']][2]
  }
  
  
  res$dependent.var = dependent.var
  return(res)
}


gems.anova.smape <- function(dat, dependent.var){ # test smape explanation power
  
  res = data.frame(metrics = c('adj.r.squared','p.value'), microtype = NA, geotype = NA, geomicrotype = NA, geomicrotype.smape = NA, geomicrotype.state.smape = NA)
  # 1. microtype alone
  f = as.formula(paste(dependent.var, 'as.factor(microtype)',sep = '~'))
  
  lm1 = lm(f, data = dat,weights = cnt_tmcs)
  res$microtype[res$metrics=='adj.r.squared'] = summary(lm1)$adj.r.squared
  
  # 2. geotype alone
  f = as.formula(paste(dependent.var, 'as.factor(geotype)',sep = '~'))
  
  lm2 = lm(f, data = dat,weights = cnt_tmcs)
  res$geotype[res$metrics=='adj.r.squared'] = summary(lm2)$adj.r.squared
  
  # 3. microtype + geotype
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype)',sep = '~'))
  
  lm3 = lm(f, data = dat,weights = cnt_tmcs)
  res$geomicrotype[res$metrics=='adj.r.squared'] = summary(lm3)$adj.r.squared
  
  # anova 3 to 1
  res$geomicrotype[res$metrics=='p.value']= anova(lm1,lm3)[['Pr(>F)']][2]
  
  
  # 4. microtype + geotype + sampe
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype) + + frc1 + frc2 + frc3 + frc4p',sep = '~'))
  
  lm4 = lm(f, data = dat,weights = cnt_tmcs)
  res$geomicrotype.smape[res$metrics=='adj.r.squared'] = summary(lm4)$adj.r.squared
  
  # anova 4 to 3
  res$geomicrotype.smape[res$metrics=='p.value']= anova(lm3,lm4)[['Pr(>F)']][2]
  
  # 5. microtype + geotype + state + smape by frcs

    f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype) + as.factor(state) + frc1 + frc2 + frc3 + frc4p',sep = '~'))
    
    lm5 = lm(f, data = dat,weights = cnt_tmcs)
    res$geomicrotype.state.smape[res$metrics=='adj.r.squared'] = summary(lm5)$adj.r.squared
    
    # anova 5to 4
    res$geomicrotype.state.smape[res$metrics=='p.value']= anova(lm4,lm5)[['Pr(>F)']][2]

  
  
  res$dependent.var = dependent.var
  return(res)
}

gems.anova.emfd <- function(dat, dependent.var){
  res = data.frame(metrics = c('adj.r.squared','p.value'), modelyear = NA, microtype = NA, geotype = NA, geomicrotype = NA, geomicrotype.modelyear = NA)
 
  # 0. modelyear alone
  f = as.formula(paste(dependent.var, 'as.factor(modelyear)',sep = '~'))
  
  lm0 = lm(f, data = dat)
  res$modelyear[res$metrics=='adj.r.squared'] = summary(lm0)$adj.r.squared
  
   # 1. microtype alone
  f = as.formula(paste(dependent.var, 'as.factor(microtype)',sep = '~'))
  
  lm1 = lm(f, data = dat)
  res$microtype[res$metrics=='adj.r.squared'] = summary(lm1)$adj.r.squared
  
  # 2. geotype alone
  f = as.formula(paste(dependent.var, 'as.factor(geotype)',sep = '~'))
  
  lm2 = lm(f, data = dat)
  res$geotype[res$metrics=='adj.r.squared'] = summary(lm2)$adj.r.squared
  
  # 3. microtype + geotype
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype)',sep = '~'))
  
  lm3 = lm(f, data = dat)
  res$geomicrotype[res$metrics=='adj.r.squared'] = summary(lm3)$adj.r.squared
  
  # anova 3 to 1
  res$geomicrotype[res$metrics=='p.value']= anova(lm1,lm3)[['Pr(>F)']][2]
  
  
  # 5. microtype + geotype + model year
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype) + as.factor(modelyear)',sep = '~'))
  
  lm4 = lm(f, data = dat)
  res$geomicrotype.modelyear[res$metrics=='adj.r.squared'] = summary(lm4)$adj.r.squared
  
  # anova 4 to 3
  res$geomicrotype.modelyear[res$metrics=='p.value']= anova(lm3,lm4)[['Pr(>F)']][2]
  
  res$dependent.var = dependent.var
  return(res)
}


gems.anova.emfd.nomodelyear <- function(dat, dependent.var){
  res = data.frame(metrics = c('adj.r.squared','p.value'), microtype = NA, geotype = NA, geomicrotype = NA)
  
  # # 0. modelyear alone
  # f = as.formula(paste(dependent.var, 'as.factor(modelyear)',sep = '~'))
  # 
  # lm0 = lm(f, data = dat)
  # res$modelyear[res$metrics=='adj.r.squared'] = summary(lm0)$adj.r.squared
  
  # 1. microtype alone
  f = as.formula(paste(dependent.var, 'as.factor(microtype)',sep = '~'))
  
  lm1 = lm(f, data = dat)
  res$microtype[res$metrics=='adj.r.squared'] = summary(lm1)$adj.r.squared
  
  # 2. geotype alone
  f = as.formula(paste(dependent.var, 'as.factor(geotype)',sep = '~'))
  
  lm2 = lm(f, data = dat)
  res$geotype[res$metrics=='adj.r.squared'] = summary(lm2)$adj.r.squared
  
  # 3. microtype + geotype
  f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype)',sep = '~'))
  
  lm3 = lm(f, data = dat)
  res$geomicrotype[res$metrics=='adj.r.squared'] = summary(lm3)$adj.r.squared
  
  # anova 3 to 1
  res$geomicrotype[res$metrics=='p.value']= anova(lm1,lm3)[['Pr(>F)']][2]
  
  
  # # 5. microtype + geotype + model year
  # f = as.formula(paste(dependent.var, 'as.factor(microtype) + as.factor(geotype) + as.factor(modelyear)',sep = '~'))
  # 
  # lm4 = lm(f, data = dat)
  # res$geomicrotype.modelyear[res$metrics=='adj.r.squared'] = summary(lm4)$adj.r.squared
  # 
  # # anova 4 to 3
  # res$geomicrotype.modelyear[res$metrics=='p.value']= anova(lm3,lm4)[['Pr(>F)']][2]
  
  res$dependent.var = dependent.var
  return(res)
}


# ## calculate cosine similarity between tracks
# cosine_sim <- function(vec1, vec2){
#   cos = cosine(vec1, vec2)
#   return(cos)
# }
# 
# # ## calculate L2 distance
# L2_distance <- function(id1, id2, data){
#   vec1 <- data.matrix(subset(data, GEOID %in% id1)[paste0("hh", 7:19)])[,-c(14)]
#   cluster.id <- subset(data, GEOID %in% id2)$sub.cluster
#   all.vectors <- data.matrix(subset(data, sub.cluster %in% cluster.id)[paste0("hh", 7:19)])[,-c(14)]
#   if (nrow(all.vectors) > 1){
#     vec2 <- colMeans(all.vectors)
#   } else {
#     vec2 <- all.vectors
#   }
#   dist <- as.numeric(dist(rbind(vec1, vec2)))
#   return(dist)
# }


