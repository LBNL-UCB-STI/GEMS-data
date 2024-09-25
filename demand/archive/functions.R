# global variables and functions used

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
      reduced <-kmeans(cluster.df, k, iter.max = 1000,nstart=1000) 
    }
    
    if(method == 'pam'){
      reduced = pam(cluster.df,k)
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
    maxdensity = max(dat[dat$flow == maxflow, 'density'])
    
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
    
    nls.fit = nls(flow ~ phi*density*(density - gridlock.density))
    
  }
  
  return(return.res)
  
}

