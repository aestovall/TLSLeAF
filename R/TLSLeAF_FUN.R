ConvertNormalToDipAndDipDir<-function(x) {
  
  colnames(x) <- c("nX","nY","nZ")
  
  # The formula using atan2() with the swapped N.x and N.y already
  # gives the correct results for facets with the normal pointing
  # upwards, so just use the sign of N.z to invert the normals if they
  # point downwards.
  
  Nsign <- as.numeric(rep(1,nrow(x)))
  Nsign[x$nZ < 0] <- -1.0
  
  dipDir_rad <- atan2(Nsign*as.numeric(x$nX), Nsign*as.numeric(x$nY))
  
  # Dip direction is measured in 360 degrees, generally clockwise from North
  dipDir_rad[dipDir_rad < 0] <- dipDir_rad[dipDir_rad < 0] + 2*pi
  
  # acos() returns values in [0, pi] but using abs() all the normals
  # are considered pointing upwards, so the actual result will be in
  # [0, pi/2] as required by the definition of dip.
  # We skip the division by r because the normal is a unit vector.
  dip_rad <- acos(abs(x$nZ))
  
  dipDir_deg = dipDir_rad * (180/pi)
  dip_deg = dip_rad * (180/pi)
  
  return(data.frame(dipDir_deg,
                    dip_deg)
  )
  
}


RZA<-function (dat, deg = FALSE){
  colnames(dat)[1:3] <- c("X","Y","Z")
  r <- sqrt(dat$X^2 + dat$Y^2 +dat$Z^2)
  inc <- ((acos(dat$Z/r)))
  remove(r)
  gc()
  
  # az <- atan2(scan$y,scan$x)
  
  if(deg==TRUE){
    inc<-inc*(180/pi)
    return(inc)
    # scan$az<-scan$az*(180/pi)
  }
  return(inc)
  
}

# RZA<-function (scan, deg = FALSE){
#   colnames(scan)[1:3] <- c("x","y","z")
#   scan$r <- sqrt(scan$x^2 + scan$y^2 +scan$z^2)
#   scan$inc <- ((acos(scan$z/scan$r)))
#   scan$az <- atan2(scan$y,scan$x)
#   
#   if(deg==TRUE){
#     scan$inc<-scan$inc*(180/pi)
#     scan$az<-scan$az*(180/pi)
#   }
#   
#   return(scan)
# }

scatter<-function(dat){
  # time<-Sys.time()
  dat[,12:14]<-dat[,1:3]/sqrt(rowSums(dat[,1:3]^2))
  dat$dot<-geometry::dot(dat[,12:14],dat[,c('nX','nY','nZ')],d=2)
  my.dt <- as.data.table(cbind(abs(dat$dot),1))
  dat$dot<-my.dt[,row.min:=pmin(V1,V2)]$row.min
  return(acos(dat$dot) * (180/pi))
  # print(Sys.time()-time)#
}

normalize_topography<-function(dat, res = 5){
  defaultW <- getOption("warn")
  options(warn = -1)
  crs <- sp::CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs")
  colnames(dat)[1:3]<-c("X","Y","Z")
  las<-LAS(dat[,1:3])
  projection(las)<-crs
  
  r <- raster(xmn=-200, xmx=200, ymn=-200, ymx=200, resolution = res)
  topo<-grid_metrics(las, quantile(Z, 0.01), r)
  # plot(topo, col = viridis::viridis(250))
  
  crs(topo)<-crs
  slope<-terrain(topo, opt = "slope", unit = "degrees", neighbors = 8)
  # plot(slope)
  
  topo[slope>40]<-NA
  setMinMax(topo)
  
  topo.df<-as.data.frame(rasterToPoints(topo))
  colnames(topo.df)<-c("X","Y","Z")
  
  ws <- seq(3,12, 3)
  th <- seq(0.1, 1.5, length.out = length(ws))
  
  topo.las<-LAS(topo.df)
  crs(topo.las)<-crs
  
  ground<-lasground(topo.las, pmf(ws, th), last_returns = FALSE)
  
  topo.las.r<-grid_terrain(ground, res = res, knnidw(k = 21))
  # plot(topo.las.r)
  
  r.p<-topo.las.r
  r.p[!is.na(r.p)]<-1
  
  topo.las.p<-rasterToPolygons(r.p, dissolve = TRUE)
  # lines(topo.las.p)
  
  lasclipPolygon(las, topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,1], topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,2])
  
  las<- las - topo.las.r
  options(warn = defaultW)
  return(las@data$Z)
}

LAvoxel<-function(dat,res = 0.1){
  defaultW <- getOption("warn")
  options(warn = -1)
  
  colnames(dat)[1:3]<-c("X","Y","Z")
  dat<-dat[dat$class==0,]
  las<-LAS(dat[,1:3])
  crs(las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  las@data$Z<-dat$z_cor
  las@data$dip_deg<-dat$dip_deg
  las@data$count<-1
  
  voxels<-voxel_metrics(las, ~mean(dip_deg), res=res)
  colnames(voxels)[4]<-"dip_dir"
  voxels$dip_dir_sd<-voxel_metrics(las, ~sd(dip_deg), res=res)[,4]
  voxels$n<-voxel_metrics(las, ~sum(count), res=res)[,4]
  
  options(warn = defaultW)
  return(voxels)
}

sim_LAD<-function(voxels){
  vox_ls<-list()
  for(i in 1:10){
    voxels.dup<-voxels
    voxels.dup$a<-apply(voxels,1,FUN = function(x)  rnorm(1, mean=x[4], sd=x[5]))
    vox_ls[[i]]<-voxels.dup
  }
  return(do.call(rbind, vox_ls))
}

# sim_LAD<-function(voxels){
#   do.call(rbind,
#           apply(voxels, 1,
#                 function(x) data.frame(X=x[1],
#                                        Y=x[2],
#                                        Z=x[3],
#                                        a=rnorm(10,
#                                                mean=x[4],
#                                                sd=x[5])))
#   )
# }

normalCalc<-function(input_file){
  
  print(paste("Input gridded point cloud and estimate normals"))
  print(paste("Processing", input_file))
  
  term<-1
  
  term<-run(paste(cloudcompare, # call Cloud Compare. The .exe file folder must be in the system PATH
                  "-SILENT",
                  "-C_EXPORT_FMT", "ASC", "-PREC", 6, #Set asc as export format
                  "-NO_TIMESTAMP",
                  "-COMPUTE_NORMALS",
                  "-O", input_file, #open the subsampled file
                  "-SAVE_CLOUDS",
                  sep = " "))
  
  while (term==1) sys.sleep(10)
}


angleCalc<-function(dat, center, scatterLim=85){
  #column naming
  colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")
  
  #Adjust coordinates to scanner center
  # dat[,1:3]<-dat[,1:3]-t(c(center))
  
  #ensures all data are numeric, avoiding errors
  dat<-as.data.frame(dat)
  for(ii in 1:length(dat)) if (class(dat[[1,ii]])=="character") dat[,ii]<- as.numeric( dat[,ii] )
  
  #calculate the radius, zenith, and azimuth angles from XYZ coordinates
  dat$inc<-RZA(dat[,1:3])
  
  # estimate scattering angle and filter, removing steep angles
  dat$scatter<-scatter(dat)
  dat<-dat[dat$scatter<=scatterLim,]
  dat<-na.omit(dat)
  
  #Convert normals to leaf orientation and leaf angle
  dat<-cbind(dat[,1:3],
             ConvertNormalToDipAndDipDir(dat[,8:10])[,2])
  
  gc()
  
  return(dat[!is.nan(dat[,4]),])
  
}

cloudRotation<-function(angle.file.name, transformationMatrix){
  print(paste("Rotating", angle.file.name))
  
  write.table(transformationMatrix, "transformation_temp.txt", 
              row.names = FALSE, col.names = FALSE)
  
  term<-1
  
  term<-run(paste(cloudcompare, # call Cloud Compare. The .exe file folder must be in the system PATH
                  "-SILENT",
                  "-C_EXPORT_FMT", "ASC", "-PREC", 6, #Set asc as export format
                  "-NO_TIMESTAMP",
                  "-AUTO_SAVE OFF",
                  "-O", angle.file.name, #open the subsampled file
                  "-APPLY_TRANS", "transformation_temp.txt",
                  "-SAVE_CLOUDS", "FILE", gsub(".asc","_t.asc",angle.file.name),
                  # "-SS SPATIAL", SS,
                  # "-OCTREE_NORMALS", scales[1],
                  # "-SAVE_CLOUDS","FILE", gsub(".asc","_0_10_NORM.asc",c2c.file),
                  # "-OCTREE_NORMALS", scales[2],
                  # "-SAVE_CLOUDS","FILE", gsub(".asc","_0_50_NORM.asc",c2c.file),
                  # "-OCTREE_NORMALS", scales[3],
                  # "-SAVE_CLOUDS", "FILE", gsub(".asc","_0_75_NORM.asc",c2c.file),
                  sep = " "))  
  
  while (term==1) sys.sleep(10)
  
}

classMetricCalc<-function(c2c.file, SS=0.02, scales=c(0.1,0.5,0.75)){
  print(paste("Processing", c2c.file))
  
  term<-1
  
  term<-run(paste(cloudcompare, # call Cloud Compare. The .exe file folder must be in the system PATH
                  "-SILENT",
                  "-C_EXPORT_FMT", "ASC", "-PREC", 6, #Set asc as export format
                  "-NO_TIMESTAMP",
                  "-AUTO_SAVE OFF",
                  "-O", c2c.file, #open the subsampled file
                  "-SS SPATIAL", SS,
                  "-OCTREE_NORMALS", scales[1],
                  "-SAVE_CLOUDS","FILE", gsub(".asc","_0_10_NORM.asc",c2c.file),
                  "-OCTREE_NORMALS", scales[2],
                  "-SAVE_CLOUDS","FILE", gsub(".asc","_0_50_NORM.asc",c2c.file),
                  "-OCTREE_NORMALS", scales[3],
                  "-SAVE_CLOUDS", "FILE", gsub(".asc","_0_75_NORM.asc",c2c.file),
                  sep = " "))  
  
  while (term==1) sys.sleep(10)
  
}

rfPrep<-function(c2c.file){
  remove(dat)
  gc()
  dat<-data.table::fread(gsub(".asc","_0_10_NORM.asc",c2c.file),header = FALSE)
  colnames(dat)[1:7]<-c("X","Y","Z","angle","nX10","nY10","nZ10")
  
  dat1<-data.table::fread(gsub(".asc","_0_50_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
  colnames(dat1)[1:3]<-c("nX50","nY50","nZ50")
  
  dat2<-data.table::fread(gsub(".asc","_0_75_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
  colnames(dat2)[1:3]<-c("nX75","nY75","nZ75")
  
  dat<-cbind(dat,dat1,dat2)
  remove(dat1,dat2)
  gc()
  
  #ensures all imports are numeric, avoiding errors
  dat<-as.data.frame(dat)
  for(ii in 1:length(dat)) if (class(dat[[1,ii]])=="character") dat[,ii]<- as.numeric( dat[,ii] )
  return(dat)
}

# TLSLeAF<-function(input_file,
#                   overwrite=TRUE,
#                   center, 
#                   scatterLim=85,
#                   SS=0.02, 
#                   scales=c(0.1,0.5,0.75),
#                   rf_model,
#                   correct.topography=TRUE,
#                   voxRes,
#                   minVoxDensity=5,
#                   superDF=FALSE,
#                   clean=TRUE,...){
#   
#   
#   #names of output files
#   output_file = gsub(".ptx",".asc", input_file)
#   angle.file.name<-gsub(".ptx","_angles.asc", input_file)
#   c2c.file<-angle.file.name
#   class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
#   # c2c.file<-gsub(".asc","_t.asc",angle.file.name)
#   gc()
#   
#   #Clear any open terminals
#   # lapply(rstudioapi::terminalList(),function(x) rstudioapi::terminalKill(x))
#   
#   #Calculate normals for gridded TLS point cloud
#   if(!file.exists(output_file)|
#      overwrite) normalCalc(input_file)
#   
#   while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
#     rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
#     Sys.sleep(10)
#   }
#   print("Terminal available")
#   # lapply(rstudioapi::terminalList(),function(x) rstudioapi::terminalKill(x))
#   
#   # if(file.size(output_file)<file.size(input_file)){
#   #   warning(paste0("Normals calculation may have been interrupted. Please check ",
#   #                  output_file)) 
#   # } 
#   
#   #Calculate scattering angle and leaf angle
#   if(!file.exists(angle.file.name)|
#      overwrite){
#     dat<-data.table::fread(output_file, header = FALSE)
#     dat<-angleCalc(dat,
#                    center, 
#                    scatterLim)
#     fwrite(dat, file = angle.file.name, sep = " ", row.names = FALSE)
#   } 
#   # else dat<-data.table::fread(angle.file.name, header = FALSE)
#   
#   # while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
#   #   rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
#   #   Sys.sleep(10)
#   # }
#   # print("Terminal available")
#   
#   
#   # if(!file.exists(gsub(".asc","_t.asc",angle.file.name))|
#   #    overwrite){
#   #   cloudRotation(angle.file.name, 
#   #                 transformationMatrix)
#   # }
#   
#   # while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
#   #   rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
#   #   Sys.sleep(10)
#   # }
#   # print("Terminal available")
#   
#   #Classify wood and leaf from random forest classifier
#   if(file.exists(c2c.file)&
#      !(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
#        file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
#        file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))|
#      overwrite) classMetricCalc(c2c.file, SS, scales)
#   
#   while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
#     rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
#     Sys.sleep(10)
#   }
#   print("Terminal available")
#   
#   if(!file.exists(class.file.name)&
#      file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
#      file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
#      file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file))&
#      (round(file.size(gsub(".asc","_0_75_NORM.asc",c2c.file)), -7)==
#      round(file.size(gsub(".asc","_0_50_NORM.asc",c2c.file)), -7))) {
#     
#       dat<-na.omit(rfPrep(c2c.file))
#       dat$predict<-predict(rf_model, dat)
#       dat<-cbind(dat[,1:4], dat$predict)
#       colnames(dat)[4:5]<-c("dip_deg", "class")
#       
#       fwrite(dat, file = class.file.name, sep = " ", row.names = FALSE)
#     
#   } else Sys.sleep(10)
#   
#   if(file.exists(class.file.name)) dat<-fread(class.file.name)
#   
#   #Correct topography and calculate the LAD and vertical LAD
#   if(correct.topography) dat$z_cor<-normalize_topography(dat) else dat$z_cor<-dat$Z
#   
#   # leaf angle voxelation and density normalization
#   voxels<-LAvoxel(na.omit(dat[dat$z_cor>1,]), voxRes)
#   
#   voxels<-voxels[voxels$n>quantile(voxels$n,0.01),]
#   voxels<-voxels[voxels$n>minVoxDensity,]
#   # voxels<-voxels[voxels$Z>1,]
#   
#   #simulate LAD from voxel statistics
#   LAD<-sim_LAD(voxels)
#   LAD<-na.exclude(LAD[LAD$a>=0 & LAD$a<=90,])
#   LAD_lim<-data.frame(a=LAD$a)
#   
#   #fit beta function and get beta parameters from LAD
#   #FIT BETA DISTRIBUTION
#   
#   m<-fitdist(as.numeric(LAD_lim$a)/90,
#              'beta',
#              method='mme')
#   
#   # Get alpha and beta parametrs
#   alpha0 <- m$estimate[1] # parameter 1
#   beta0 <- m$estimate[2] # parameter 2
#   beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01),
#                                                    alpha0,beta0))
#   param<-data.frame(alpha0,beta0)
#   
#   TLSLeAF.dat<-new("TLSLeAF",
#                    parameters=data.frame(c(file=input_file,
#                                            center,
#                                            scatterLim=85,
#                                            SS=0.02,
#                                            scales=c(0.1,0.5,0.75),
#                                            voxRes=voxRes,
#                                            superDF=TRUE)),
#                    dat=dat,
#                    voxels=as.data.frame(voxels),
#                    LAD=LAD,
#                    Beta_parameters=param,
#                    beta=beta,
#                    G=data.frame(inc_bin=1, G_p=1,assumption="",density=""))
#   
#   TLSLeAF.dat@G<-G_calculations(TLSLeAF.dat)
#   
#   if(superDF) return(TLSLeAF.dat)
#   
#   if(clean==TRUE){
#     clean.temp(output_file,c2c.file)
#   }
# }
TLSLeAF<-function(input_file,
                  overwrite=TRUE,
                  center, 
                  scatterLim=85,
                  SS=0.02, 
                  scales=c(0.1,0.5,0.75),
                  # rf_model,
                  rf_model_path,
                  correct.topography=TRUE,
                  voxRes,
                  minVoxDensity=5,
                  superDF=FALSE,
                  clean=TRUE,...){
  
  
  #names of output files
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  # c2c.file<-gsub(".asc","_t.asc",angle.file.name)
  gc()
  
  #Clear any open terminals
  # lapply(rstudioapi::terminalList(),function(x) rstudioapi::terminalKill(x))
  
  #Calculate normals for gridded TLS point cloud
  if(!file.exists(output_file)|
     overwrite) normalCalc(input_file)
  
  while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
    rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
    Sys.sleep(10)
  }
  # print("Terminal available")
  # lapply(rstudioapi::terminalList(),function(x) rstudioapi::terminalKill(x))
  
  # if(file.size(output_file)<file.size(input_file)){
  #   warning(paste0("Normals calculation may have been interrupted. Please check ",
  #                  output_file)) 
  # } 
  
    #Calculate scattering angle and leaf angle
    if((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
      rstudioapi::terminalBusy(rstudioapi::terminalList())])==0)&
      (!file.exists(angle.file.name)|
       overwrite)){
      remove(dat)
      gc()
      
      dat<-data.table::fread(output_file, header = FALSE, 
                             colClasses = rep('numeric',10))
      dat<-angleCalc(dat,
                     center,
                     scatterLim)
      
      fwrite(dat, file = angle.file.name, sep = " ", row.names = FALSE)
      
      
      # Sys.sleep(20)
    }
  # else dat<-data.table::fread(angle.file.name, header = FALSE)
  
  # while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
  #   rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
  #   Sys.sleep(10)
  # }
  # print("Terminal available")
  
  
  # if(!file.exists(gsub(".asc","_t.asc",angle.file.name))|
  #    overwrite){
  #   cloudRotation(angle.file.name, 
  #                 transformationMatrix)
  # }
  
  # while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
  #   rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
  #   Sys.sleep(10)
  # }
  # print("Terminal available")
  
  #Classify wood and leaf from random forest classifier
  if(file.exists(angle.file.name)&
     !(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))|
     overwrite) classMetricCalc(angle.file.name, SS, scales)
  
  while (length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
    rstudioapi::terminalBusy(rstudioapi::terminalList())])>0) {
    Sys.sleep(10)
  }
  print("Terminal available")
  
  if(!file.exists(class.file.name)&
     file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
     file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
     file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file))&
     (round(file.size(gsub(".asc","_0_75_NORM.asc",c2c.file)), -7)==
      round(file.size(gsub(".asc","_0_50_NORM.asc",c2c.file)), -7))) {
    
    
    remove(dat)
    gc()
    
    #load the randomForest model
    dat<-na.omit(rfPrep(c2c.file))
    
    rf_model<-readRDS(rf_model_path)
    dat$predict<-predict(rf_model, dat)
    remove(rf_model)
    gc()
    
    dat<-cbind(dat[,1:4], dat$predict)
    colnames(dat)[4:5]<-c("dip_deg", "class")
    gc()
    
    fwrite(dat, file = class.file.name, sep = " ", row.names = FALSE)
    Sys.sleep(20)
    
  } else if(file.exists(class.file.name)) dat<-fread(class.file.name)
  
  #Correct topography and calculate the LAD and vertical LAD
  if(correct.topography) dat$z_cor<-normalize_topography(dat) else dat$z_cor<-dat$Z
  
  # if(colnames(dat)[6]=="z_cor"){
  
  # leaf angle voxelation and density normalization
  voxels<-LAvoxel(na.omit(dat[dat$z_cor>1,]), voxRes)
  
  voxels<-voxels[voxels$n>quantile(voxels$n,0.01),]
  voxels<-voxels[voxels$n>minVoxDensity,]
  # voxels<-voxels[voxels$Z>1,]
  
  #simulate LAD from voxel statistics
  LAD<-sim_LAD(voxels)
  LAD<-na.exclude(LAD[LAD$a>=0 & LAD$a<=90,])
  LAD_lim<-data.frame(a=LAD$a)
  
  #fit beta function and get beta parameters from LAD
  #FIT BETA DISTRIBUTION
  
  m<-fitdist(as.numeric(LAD_lim$a)/90,
             'beta',
             method='mme')
  
  # Get alpha and beta parametrs
  alpha0 <- m$estimate[1] # parameter 1
  beta0 <- m$estimate[2] # parameter 2
  beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01),
                                                   alpha0,beta0))
  param<-data.frame(alpha0,beta0)
  
  remove(m)
  gc()
  
  TLSLeAF.dat<-new("TLSLeAF",
                   parameters=data.frame(c(file=input_file,
                                           center,
                                           scatterLim=85,
                                           SS=0.02,
                                           scales=c(0.1,0.5,0.75),
                                           voxRes=voxRes,
                                           superDF=TRUE)),
                   dat=dat,
                   voxels=as.data.frame(voxels),
                   LAD=LAD,
                   Beta_parameters=param,
                   beta=beta,
                   G=data.frame(inc_bin=1, G_p=1,assumption="",density=""))
  
  remove(dat)
  gc()
  
  G<-G_calculations(TLSLeAF.dat)
  TLSLeAF.dat@G<-G
  
  if(superDF) return(TLSLeAF.dat)
  
  if(clean==TRUE){
    clean.temp(output_file,c2c.file)
  }
  # } else sys
}

TLSLeAF.class<-setClass("TLSLeAF",representation=representation(
  parameters = "data.frame",
  dat = "data.frame",
  voxels="data.frame",
  LAD = 'data.frame',
  Beta_parameters = "data.frame",
  beta="data.frame",
  G="data.frame"
))

clean.temp<-function(output_file,c2c.file){
  temp.files<-c(output_file, gsub(".asc","_angles.asc",output_file),
                gsub(".asc","_0_10_NORM.asc",c2c.file),
                gsub(".asc","_0_50_NORM.asc",c2c.file),
                gsub(".asc","_0_75_NORM.asc",c2c.file))[file.exists(c(output_file, gsub(".asc","_angles.asc",output_file),
                                                                      gsub(".asc","_0_10_NORM.asc",c2c.file),
                                                                      gsub(".asc","_0_50_NORM.asc",c2c.file),
                                                                      gsub(".asc","_0_75_NORM.asc",c2c.file)))]
  
  if(clean & length(temp.files)>0) file.remove(temp.files)
  
  gc()
}


pgapF<-function (scan.sub, row_res, col_res, z_res) {
  
  total<-(row_res/150)*col_res*5
  list <- seq(floor(min(scan.sub[3])), 
              ceiling(max(scan.sub[3])), 
              by = z_res)
  
  pgap <- NULL
  for (i in 1:length(list)) {
    pgap[i] <- 1-(length(subset(scan.sub$z, scan.sub$z<list[i])))/total
  }
  pgap<- (1-max(pgap))+pgap
  return(pgap)
}

a_vai<-function(scan, a_ls, z_res, col_res, row_res){
  vai_a<-list()
  for (a in 1:length(a_ls)){
    
    scan.sub<-subset(scan,scan$inc>a_ls[a]&scan$inc<(a_ls[a]+5))
    list <- seq(floor(min(scan.sub[3])), 
                ceiling(max(scan.sub[3])), 
                by = z_res)
    
    pgap<-pgapF(scan.sub, row_res, col_res, z_res)
    
    # pgap[pgap<exp(8/-1.1)] <- exp(8/-1.1)
    
    f <- -1.1*log(pgap)
    # 0.4/cos(a_ls[a]*(pi/180))
    
    spl <- smooth.spline(x=list,y=f, df=length(list)-round(0.3*length(list)))
    plot(spl)
    pred <- predict(spl)
    
    f.prime <- (diff(f)/diff(list))
    f.pred.prime <- predict(spl, deriv=1)
    
    f.pred.prime$y[f.pred.prime$y<0]<-0
    plot(f.pred.prime)
    lines(x = f.pred.prime, col = "red")
    
    f.prime.df<-data.frame(y=f.pred.prime$y, x=f.pred.prime$x)
    f.prime.sub<-subset(f.prime.df,f.prime.df$x<100)
    f.prime.df$y[f.pred.prime$y<0]<-0
    
    vai_a[[a]]<-data.frame(
      # scan = gsub(".ptx","", files[l]),
      a = a_ls[a],
      z = f.prime.sub[,2],
      dx = f.prime.sub[,1])
    
    
  }
  return(data.frame(do.call(rbind,vai_a), scan = name))
}

vai<-function (scan, a_range, a_bin, z_res, col_res, row_res, name){
  # scan<-RZA(scan)
  a_ls<-seq(min(a_range),max(a_range), by = a_bin) # zenith angle range
  
  vai_out<-a_vai(scan, a_ls, z_res, col_res, row_res)
  return(vai_out)
}



vai_weighted <- function(a_vai, adj_z, z_res) {
  a_vai <- a_vai[a_vai$z>0,]
  a_vai$z_adj<-a_vai$z+adj_z
  a_vai$z_bin<-floor(a_vai$z_adj)
  #weight the profiles by zenith angle
  a_vai$weight<-sin(a_vai$a*(pi/180))
  a_vai_weight<-a_vai %>%
    group_by(scan, z_bin) %>% 
    mutate(dxw = weighted.mean(dx, weight))
  return(a_vai_weight)
}

G_calculations<-function(df){
  
  G_ls<-list()
  dat<-df@dat
  dat<-dat[dat$class==0,]
  dat$inc<-RZA(dat[,1:3], deg = TRUE)
  dat$inc_bin<-floor(dat$inc)
  dat$dip_bin<-floor(dat$dip_deg)
  dat<-dat[dat$inc_bin>=0&dat$inc_bin<=90,]
  
  dat$count<-1
  
  dat.sub<-rbind(data.frame(inc_bin=dat$inc_bin,
                            dip_bin=dat$dip_bin,
                            count=dat$count),
                 data.frame(inc_bin=0:90,
                            dip_bin=0:90,
                            count=0))
  
  G_test<-aggregate(count~dip_bin, FUN="sum",dat.sub)
  G_test$p<-G_test$count/sum(G_test$count)
  plot(G_test)
  
  ran_G_ls<-list()
  for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(G_test[G_test$dip_bin==c(0:90)[i],],
                                                     inc_bin=0:90, row.names = NULL)
  G_test<-do.call(rbind, ran_G_ls)
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"random"
  G_calc$density<-"non-normalized"
  
  plot(G_calc$inc_bin,
       G_calc$G_p)
  
  G_ls[[1]]<-G_calc
  
  G_test<-aggregate(count~inc_bin+dip_bin, FUN="sum",dat.sub)
  p<-aggregate(count~inc_bin, FUN="sum",
               G_test)
  
  G_test<-merge(G_test,p, by="inc_bin")
  
  G_test$p<-G_test$count.x/G_test$count.y
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"non-random"
  G_calc$density<-"non-normalized"
  
  G_ls[[2]]<-G_calc
  
  dat<-df@voxels
  dat$inc<-RZA(df@voxels[,1:3], deg = TRUE)
  dat$inc_bin<-floor(dat$inc)
  dat$dip_bin<-floor(dat$dip_dir)
  dat$count<-1
  
  dat.sub<-rbind(data.frame(inc_bin=dat$inc_bin,
                            dip_bin=dat$dip_bin,
                            count=dat$count),
                 data.frame(inc_bin=0:90,
                            dip_bin=0:90,
                            count=0))
  
  G_test<-aggregate(count~dip_bin, FUN="sum",dat.sub)
  G_test$p<-G_test$count/sum(G_test$count)
  # plot(G_test)
  
  ran_G_ls<-list()
  for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(G_test[G_test$dip_bin==c(0:90)[i],],
                                                     inc_bin=0:90, row.names = NULL)
  G_test<-do.call(rbind, ran_G_ls)
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"random"
  G_calc$density<-"normalized"
  
  G_ls[[3]]<-G_calc
  
  G_test<-aggregate(count~inc_bin+dip_bin, FUN="sum",dat.sub)
  p<-aggregate(count~inc_bin, FUN="sum",
               G_test)
  
  G_test<-merge(G_test,p, by="inc_bin")
  
  G_test$p<-G_test$count.x/G_test$count.y
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"non-random"
  G_calc$density<-"normalized"
  
  G_ls[[4]]<-G_calc
  
  return(do.call(rbind, G_ls))
}

ptx.header<-function(input_file){
  list(
    name = scan,
    col_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, header = FALSE)[1,]),
    row_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, header = FALSE)[2,]),
    scan_center = as.matrix(read.csv(input_file, sep = "", skip = 2,nrows = 1, header = FALSE)),
    reg = as.matrix(read.csv(input_file, sep = "", skip = 3,nrows = 3, header = FALSE)),
    trans_matrix = as.matrix(read.csv(input_file, sep = "", skip = 6,nrows = 4, header = FALSE))
  )
}

pgap.angle<-function(scan,
                     z_res=1,
                     adj_z=1.6,
                     a_range=c(10,60),
                     a_bin = 5,
                     plane = FALSE,
                     header,...){
  
  
  
  #calculate inclination angle
  inc<-RZA(scan[,1:3], deg=TRUE)
  scan$inc<-inc
  remove(inc)
  
  # zenith angle range
  a_ls<-seq(min(a_range),max(a_range), by = a_bin) 
  
  pgap_ls<-list()
  
  for (a in 1:length(a_ls)){
    # test_ls<-list()
    
    scan.sub<-subset(scan,scan$inc>a_ls[a]&scan$inc<(a_ls[a]+5))
    scan.sub<-na.exclude(scan.sub)
    # if(nrow(scan.sub) ==0) next
    list <- seq(floor(min(scan.sub[3])), 
                ceiling(max(scan.sub[3])), 
                by = z_res)
    
    pgap<-pgapF(scan.sub, header$row_res, header$col_res, z_res)
    # plot(pgap)
    
    pgap_ls[[a]]<-data.frame(list,pgap, inc_bin=a_ls[a])
  }
  remove(scan)
  remove(scan.sub)
  gc()
  pgap_all<-do.call(rbind, pgap_ls)
  return(pgap_all)
}


pgap2PAI <- function(pgap_all, method=NULL, G=NULL,...){
  
  if((method != "A"| 
      method != "D") & !is.null(G)) G<-G[G$inc_bin %in% unique(pgap_all$inc_bin),]
  
  if( (method != "A"| 
       method != "D") & is.null(G) ) print("No G provided!")
  
  if(method=="A"){
    pgap_all$f <- cos(pgap_all$inc_bin*(pi/180))/0.5*(-1)*log(pgap_all$pgap)
    pgap_all$G<-0.5
  } 
  
  if(method=="B"){
    pgap_all<-merge(pgap_all, G, by.x="inc_bin")
    pgap_all$f <- cos(pgap_all$inc_bin*(pi/180))/
      pgap_all$G*(-1)*log(pgap_all$pgap)
  } 
  if(method=="C"){
    pgap_all<-merge(pgap_all, G, by.x="inc_bin")
    pgap_all$f <- cos(pgap_all$inc_bin*(pi/180))/
      pgap_all$G*(-1)*log(pgap_all$pgap)
  } 
  
  if(method=="D"){
    pgap_all$f <- -1*log(pgap_all$pgap)
    pgap_all$x_theta<-2/pi*tan(pgap_all$inc_bin*(pi/180))
    
    temp.ls<-list()
    for(j in 1:length(unique(pgap_all$list))){
      temp<-pgap_all[pgap_all$list==unique(pgap_all$list)[j],]
      temp.ls[[j]]<-data.frame(list = unique(list)[j],
                               f=sum(coef(lm(f~x_theta,temp))),
                               sd=sigma(lm(f~x_theta,temp)),
                               se=sigma(lm(f~x_theta,temp))/sqrt(nrow(temp)))
    }
    pgap_lm<-do.call(rbind, temp.ls)
    pgap_lm$f.cimin<-pgap_lm$f-1.96*pgap_lm$se
    pgap_lm$f.cimax<-pgap_lm$f+1.96*pgap_lm$se
    
    pgap_lm<-pgap_lm[order(pgap_lm$list),]
    pgap_lm$sd[pgap_lm$list> na.exclude(pgap_lm$list[pgap_lm$f==max(pgap_lm$f, na.rm=TRUE)])]<-0
    pgap_lm$f.cimin[pgap_lm$list> na.exclude(pgap_lm$list[pgap_lm$f==max(pgap_lm$f, na.rm=TRUE)])]<-0
    pgap_lm$f.cimax[pgap_lm$list> na.exclude(pgap_lm$list[pgap_lm$f==max(pgap_lm$f, na.rm=TRUE)])]<-0
    
    pgap_lm$f[pgap_lm$list> na.exclude(pgap_lm$list[pgap_lm$f==max(pgap_lm$f, na.rm=TRUE)])]<-max(pgap_lm$f, na.rm=TRUE)
    
    pgap_lm$sd[is.na(pgap_lm$f)]<-0
    pgap_lm$f.cimin[is.na(pgap_lm$f)]<-0
    pgap_lm$f.cimax[is.na(pgap_lm$f)]<-0
    pgap_lm$f[is.na(pgap_lm$f)]<-0
    
    pgap_all<-pgap_lm
    
    pgap_all$f.prime.ci_min<-c((diff(pgap_all$f.cimin)/diff(pgap_all$list)),0)
    pgap_all$f.prime.ci_max<-c((diff(pgap_all$f.cimax)/diff(pgap_all$list)),0)
    pgap_all$f.prime.ci_min[pgap_all$f.prime==0]<-0
    pgap_all$f.prime.ci_max[pgap_all$f.prime==0]<-0
    
    pgap_all$f.prime<-c((diff(pgap_all$f)/diff(pgap_all$list)),0)
    
    spl <- smooth.spline(x=pgap_all$list,
                         y=pgap_all$f, 
                         df=length(pgap_all$list)-round(0.3*length(pgap_all$list)))
    pgap_all$pred <- predict(spl)$y
    pgap_all$f.pred.prime <- predict(spl, deriv=1)$y
    pgap_all$f.pred.prime[pgap_all$f.pred.prime<0]<-0
    return(pgap_all)
  } 
  
  if(method!="D"){
    pgap_all<-lapply(unique(pgap_all$inc_bin), function(x){
      pgap.sub<-pgap_all[pgap_all$inc_bin==x,]
      pgap.sub$f.prime<-c((diff(pgap.sub$f)/diff(pgap.sub$list)),0)
      
      spl <- smooth.spline(x=pgap.sub$list,
                           y=pgap.sub$f, 
                           df=length(pgap.sub$list)-round(0.3*length(pgap.sub$list)))
      pgap.sub$pred <- predict(spl)$y
      pgap.sub$f.pred.prime <- predict(spl, deriv=1)$y
      pgap.sub$f.pred.prime[pgap.sub$f.pred.prime<0]<-0
      return(pgap.sub)
    })
    
    return(do.call(rbind,pgap_all)) }
}
  

G_calculations<-function(df){
  
  G_ls<-list()
  dat<-df@dat
  dat<-dat[dat$class==0,]
  dat$inc<-RZA(dat[,1:3], deg = TRUE)
  dat$inc_bin<-floor(dat$inc)
  dat$dip_bin<-floor(dat$dip_deg)
  dat<-dat[dat$inc_bin>=0&dat$inc_bin<=90,]
  
  dat$count<-1
  
  dat.sub<-rbind(data.frame(inc_bin=dat$inc_bin,
                            dip_bin=dat$dip_bin,
                            count=dat$count),
                 data.frame(inc_bin=0:90,
                            dip_bin=0:90,
                            count=0))
  
  G_test<-aggregate(count~dip_bin, FUN="sum",dat.sub)
  G_test$p<-G_test$count/sum(G_test$count)
  # plot(G_test)
  
  ran_G_ls<-list()
  for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(G_test[G_test$dip_bin==c(0:90)[i],],
                                                     inc_bin=0:90)
  G_test<-do.call(rbind, ran_G_ls)
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"random"
  G_calc$density<-"non-normalized"
  
  plot(G_calc$inc_bin,
       G_calc$G_p)
  
  G_ls[[1]]<-G_calc
  
  G_test<-aggregate(count~inc_bin+dip_bin, FUN="sum",dat.sub)
  p<-aggregate(count~inc_bin, FUN="sum",
               G_test)
  
  G_test<-merge(G_test,p, by="inc_bin")
  
  G_test$p<-G_test$count.x/G_test$count.y
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"non-random"
  G_calc$density<-"non-normalized"
  
  G_ls[[2]]<-G_calc
  
  dat<-df@voxels
  dat$inc<-RZA(df@voxels[,1:3], deg = TRUE)
  dat$inc_bin<-floor(dat$inc)
  dat$dip_bin<-floor(dat$dip_dir)
  dat$count<-1
  
  dat.sub<-rbind(data.frame(inc_bin=dat$inc_bin,
                            dip_bin=dat$dip_bin,
                            count=dat$count),
                 data.frame(inc_bin=0:90,
                            dip_bin=0:90,
                            count=0))
  
  G_test<-aggregate(count~dip_bin, FUN="sum",dat.sub)
  G_test$p<-G_test$count/sum(G_test$count)
  # plot(G_test)
  
  ran_G_ls<-list()
  for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(G_test[G_test$dip_bin==c(0:90)[i],],
                                                     inc_bin=0:90)
  G_test<-do.call(rbind, ran_G_ls)
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"random"
  G_calc$density<-"normalized"
  
  G_ls[[3]]<-G_calc
  
  G_test<-aggregate(count~inc_bin+dip_bin, FUN="sum",dat.sub)
  p<-aggregate(count~inc_bin, FUN="sum",
               G_test)
  
  G_test<-merge(G_test,p, by="inc_bin")
  
  G_test$p<-G_test$count.x/G_test$count.y
  
  G_test$A_test<-abs(
    ( 1/tan(G_test$inc_bin*(pi/180)) )*
      ( 1/tan(G_test$dip_bin*(pi/180)) )
  )
  
  G_test$A1<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))
  
  G_test$phi<-NA
  G_test$phi[G_test$A_test<=1]<-acos(
    (1/tan(G_test$inc_bin[G_test$A_test<=1]*(pi/180)))*
      (1/tan(G_test$dip_bin[G_test$A_test<=1]*(pi/180)))
  )
  
  G_test$A2<-
    cos(G_test$inc_bin*(pi/180))*
    cos(G_test$dip_bin*(pi/180))*
    (1+(2/pi)*(tan(G_test$phi)-G_test$phi))
  
  G_test$A<-G_test$A1
  G_test$A[G_test$A_test<=1]<-G_test$A2[G_test$A_test<=1]
  
  G_test$G_p<-G_test$A*G_test$p
  
  G_calc<-aggregate(G_p~inc_bin, 
                    FUN= "sum", 
                    G_test, na.rm=TRUE)
  
  G_calc$assumption<-"non-random"
  G_calc$density<-"normalized"
  
  G_ls[[4]]<-G_calc
  
  return(do.call(rbind, G_ls))
}