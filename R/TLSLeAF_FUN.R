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


RZA<-function (scan, deg = FALSE){
  colnames(scan)[1:3] <- c("x","y","z")
  scan$r <- sqrt(scan$x^2 + scan$y^2 +scan$z^2)
  scan$inc <- ((acos(scan$z/scan$r)))
  scan$az <- atan2(scan$y,scan$x)
  
  if(deg==TRUE){
    scan$inc<-scan$inc*(180/pi)
    scan$az<-scan$az*(180/pi)
  }
  
  return(scan)
}

scatter<-function(dat){
  # time<-Sys.time()
  dat[,14:16]<-dat[,1:3]/sqrt(rowSums(dat[,1:3]^2))
  dat$dot<-geometry::dot(dat[,14:16],dat[,c('nX','nY','nZ')],d=2)
  my.dt <- as.data.table(cbind(abs(dat$dot),1))
  dat$dot<-my.dt[,row.min:=pmin(V1,V2)]$row.min
  return(acos(dat$dot) * (180/pi))
  # print(Sys.time()-time)#
}

normalize_topography<-function(dat, res = 5){
  las<-LAS(dat[,1:3])
  crs(las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  
  r <- raster(xmn=-200, xmx=200, ymn=-200, ymx=200, resolution = res)
  topo<-grid_metrics(las, quantile(Z, 0.01), r)
  # plot(topo, col = viridis::viridis(250))
  
  crs(topo)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  slope<-terrain(topo, opt = "slope", unit = "degrees", neighbors = 8)
  # plot(slope)
  
  topo[slope>40]<-NA
  setMinMax(topo)
  
  topo.df<-as.data.frame(rasterToPoints(topo))
  colnames(topo.df)<-c("X","Y","Z")
  
  ws <- seq(3,12, 3)
  th <- seq(0.1, 1.5, length.out = length(ws))
  
  topo.las<-LAS(topo.df)
  crs(topo.las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  
  ground<-lasground(topo.las, pmf(ws, th), last_returns = FALSE)
  
  topo.las.r<-grid_terrain(ground, res = res, knnidw(k = 21))
  # plot(topo.las.r)
  
  r.p<-topo.las.r
  r.p[!is.na(r.p)]<-1
  
  topo.las.p<-rasterToPolygons(r.p, dissolve = TRUE)
  # lines(topo.las.p)
  
  lasclipPolygon(las, topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,1], topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,2])
  
  las<- las - topo.las.r
  
  return(las@data$Z)
}

LAvoxel<-function(dat,res = 0.1){
  las<-LAS(dat[,1:3])
  crs(las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  las@data$Z<-dat$z_cor
  las@data$dip_deg<-dat$dip_deg
  las@data$count<-1
  
  voxels<-voxel_metrics(las, mean(dip_deg), res=res)
  colnames(voxels)[4]<-"dip_dir"
  voxels$dip_dir_sd<-voxel_metrics(las, sd(dip_deg), res=res)[,4]
  voxels$n<-voxel_metrics(las, sum(count), res=res)[,4]
  
  return(voxels)
}

sim_LAD<-function(x,sd){
  sim_ls<-list()
  list_length<-length(x)
  pb <- txtProgressBar(min = 0, max = list_length, style = 3)
  for(i in 1:list_length){
    sim_ls[[i]]<-data.frame(a=rnorm(10, mean=x[i], sd=sd[i]))
    setTxtProgressBar(pb, i)
  }
  return(data.frame(a=do.call(rbind, sim_ls)))
  close(pb)
}

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


angleCalc<-function(dat, center, SCATTER_LIM=85){
  #column naming
  colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")
  
  #Adjust coordinates to scanner center
  dat[,1:3]<-dat[,1:3]-t(c(center))
  
  #ensures all data are numeric, avoiding errors
  dat<-as.data.frame(dat)
  for(ii in 1:length(dat)) if (class(dat[[1,ii]])=="character") dat[,ii]<- as.numeric( dat[,ii] )
  
  #calculate the radius, zenith, and azimuth angles from XYZ coordinates
  dat<-RZA(dat)
  
  # estimate scattering angle and filter, removing steep angles
  dat$scatter<-scatter(dat)
  dat<-dat[dat$scatter<=SCATTER_LIM,]
  dat<-na.omit(dat)
  
  #Convert normals to leaf orientation and leaf angle
  dat<-cbind(dat[,1:3],
             ConvertNormalToDipAndDipDir(dat[,8:10])[,2])
  return(dat)

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



TLSLeAF<-function(input_file,
                  overwrite=TRUE,
                  center, 
                  SCATTER_LIM=85,
                  SS=0.02, 
                  scales=c(0.1,0.5,0.75),
                  rf_model,
                  vox.res,
                  minVoxDensity=5,
                  superDF=FALSE,
                  clean=TRUE,...){
  
  
  #names of output files
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  gc()
  
  #Calculate normals for gridded TLS point cloud
  if(!file.exists(output_file)|
     overwrite) normalCalc(input_file)
  
  #Calculate scattering angle and leaf angle
  if(!file.exists(angle.file.name)|
     overwrite){
    dat<-data.table::fread(output_file, header = FALSE)
    dat<-angleCalc(dat,
                   center, 
                   SCATTER_LIM)
    fwrite(dat, file = angle.file.name, sep = " ", row.names = FALSE)
  } 
  # else dat<-data.table::fread(angle.file.name, header = FALSE)
  
  #Classify wood and leaf from random forest classifier
  if(file.exists(c2c.file)&
      !(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))|
     overwrite) classMetricCalc(c2c.file, SS, scales)
  
  if(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file))) {
    
    dat<-na.omit(rfPrep(c2c.file))
    dat$predict<-predict(rf_model, dat)
    dat<-cbind(dat[,1:4], dat$predict)
    colnames(dat)[4:5]<-c("dip_deg", "class")
    
  }
  
  #Correct topography and calculate the LAD and vertical LAD
  if(correct.topography == TRUE) dat$z_cor<-normalize_topography(dat) else dat$z_cor<-dat$Z
  
  # leaf angle voxelation and density normalization
  voxels<-LAvoxel(dat, vox.res)
  voxels<-voxels[voxels$n>minVoxDensity,]
  
  
  #simulate LAD from voxel statistics
  LAD<-sim_LAD(voxels$dip_dir,voxels$dip_dir_sd)
  LAD<-na.exclude(LAD[LAD$a>=0 & LAD$a<=90,])
  LAD_lim<-data.frame(a=LAD)
  
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
  
  TLSLeAF.dat<-new("TLSLeAF",
                   parameters=data.frame(c(file=input_file, 
                                           center, 
                                           SCATTER_LIM=85,
                                           SS=0.02, 
                                           scales=c(0.1,0.5,0.75),
                                           vox.res=vox.res,
                                           superDF=FALSE)),
                   dat=dat,
                   voxels=as.data.frame(voxels),
                   LAD=LAD_lim,
                   Beta_parameters=param)
  
  if(superDF) return(TLSLeAF.dat)
  
}


TLSLeAF.class<-setClass("TLSLeAF",representation=representation(
  parameters = "data.frame",
  dat = "data.frame",
  voxels="data.frame",
  LAD = 'data.frame',
  Beta_parameters = "data.frame"
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


