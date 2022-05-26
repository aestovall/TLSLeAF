
ConvertNormalToDipAndDipDir<-function(x) {
  
  Nsign<-dipDir_rad<-dip_rad<-dip_deg<-NULL
  
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
  r<-inc<-NULL
  
  colnames(dat)[1:3] <- c("X","Y","Z")
  r <- sqrt(dat$X^2 + dat$Y^2 +dat$Z^2)
  
  inc<-matrix(NA, nrow=length(r))
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

scatter<-function(dat, cols_num){
  my.dt<-NULL
  # time<-Sys.time()
  # dat<-cbind(dat, )
  dat.add<-dat[,1:3]/sqrt(rowSums(dat[,1:3]^2))
  colnames(dat.add)<-c("sX","sY","sZ")
  dat<-cbind(dat,dat.add)
  remove(dat.add)
  
  dat$dot<-geometry::dot(dat[,cols_num+1:3],dat[,c('nX','nY','nZ')],d=2)
  my.dt <- as.data.table(cbind(abs(dat$dot),1))
  dat$dot<-my.dt[,row.min:=pmin(V1,V2)]$row.min
  return(acos(dat$dot) * (180/pi))
  remove(my.dt)
  # print(Sys.time()-time)#
}

normalize_topography<-function(x, res = 5){
  las<-r<-topo<-topo.df<-ground<-topo.las.r<-r.p<-topo.p<-slope<-poly.id<-NULL
  
  # defaultW <- getOption("warn")
  # options(warn = -1)
  # crs_tls <- sp::CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs")
  crs_tls<-sp::CRS("+init=epsg:26918")
  # crs_tls
  colnames(x)[1:3]<-c("X","Y","Z")
  las<-LAS(x[,1:3], proj4string = crs_tls)
  # projection(las)<-crs
  las@data$dip_deg<-x[,4]
  las@data$class<-x[,5]
  
  r <- raster(xmn=-200, xmx=200, ymn=-200, ymx=200, resolution = res)
  topo<-grid_metrics(las, quantile(Z, 0.01), r)
  crs(topo)<-crs_tls
  # plot(topo)
  
  topo.df<-as.data.frame(rasterToPoints(topo))
  colnames(topo.df)<-c("X","Y","Z")
  
  ws <- seq(3,12, 3)
  th <- seq(0.1, 1.5, length.out = length(ws))
  
  topo.las<-LAS(topo.df)
  crs(topo.las)<-crs_tls
  
  ground<-classify_ground(topo.las, pmf(ws, th), last_returns = FALSE)
  topo.las.r<-grid_terrain(ground, r, kriging(k = 10))
  
  # topo<-focal(topo, matrix(1,3,3))
  # topo.smooth<-disaggregate(topo, 5, method='bilinear')
  # # density.r<-grid_density(las, topo.smooth)
  # # plot(topo, col = viridis::viridis(250))
  # plot(topo.las.r)
  
  r.p<-buffer(topo.las.r, width=1)
  r.p[!is.na(r.p)]<-1
  
  topo.las.p<-rasterToPolygons(r.p, dissolve = TRUE)
  # lines(topo.las.p)
  
  las<-clip_polygon(las, topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,1], topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,2])
  
  # plot(las)
  las<- las - topo.las.r
  # plot(las)
  
  slope<-terrain(topo, opt = "slope", unit = "degrees", neighbors = 8)
  # plot(slope)
  
  topo[is.na(slope)]<-NA
  topo.p<-topo
  topo.p[!is.na(topo.p)]<-1
  
  topo.las.p<-rasterToPolygons(topo.p, dissolve = TRUE)
  # lines(topo.las.p)
  
  poly.id<-length(topo.las.p@polygons[[1]]@Polygons)
  
  topo.las.p@polygons[[1]]@Polygons
  
  a.ls<-c(do.call(rbind, lapply(1:poly.id, function(x){
    area.sub<-topo.las.p@polygons[[1]]@Polygons[[x]]
    return(area.sub@area)
  })))
  
  max.poly.id<-(1:poly.id)[max(a.ls)==a.ls]
  
  las<-clip_polygon(las, topo.las.p@polygons[[1]]@Polygons[[max.poly.id]]@coords[,1], topo.las.p@polygons[[1]]@Polygons[[max.poly.id]]@coords[,2])
  # plot(las)
  return(las@data)
  
  # remove(las)
  
}

LAvoxel<-function(x,res = 0.1){
  # defaultW <- getOption("warn")
  # options(warn = -1)
  las<-voxels<-NULL
  colnames(x)[1:3]<-c("X","Y","Z")
  x<-x[x$class==0,]
  crs_tls<-sp::CRS("+init=epsg:26918")
  las<-LAS(x[,1:3],proj4string = crs_tls)
  
  # crs(las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  # las@data$Z<-x$z_cor
  las@data$dip_deg<-x$dip_deg
  las@data$count<-1
  
  voxels<-voxel_metrics(las, ~mean(dip_deg), res=res)
  colnames(voxels)[4]<-"dip_dir"
  voxels$dip_dir_sd<-voxel_metrics(las, ~sd(dip_deg), res=res)[,4]
  voxels$n<-voxel_metrics(las, ~sum(count), res=res)[,4]
  
  # options(warn = defaultW)
  return(voxels)
  
  remove(las)
  remove(voxels)
}

sim_LAD<-function(x){
  voxels.dup<-NULL
  vox_ls<-list()
  for(i in 1:10){
    voxels.dup<-x
    voxels.dup$a<-apply(voxels.dup,1,FUN = function(xx)  rnorm(1, mean=xx[4], sd=xx[5]))
    vox_ls[[i]]<-voxels.dup
    remove(voxels.dup)
  }
  return(do.call(rbind, vox_ls))
  remove(vox_ls)
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
                  "-AUTO_SAVE OFF",
                  "-COMPUTE_NORMALS",
                  "-O", input_file, #open the subsampled file
                  # "-NORMALS_TO_DIP",
                  "-SAVE_CLOUDS", 
                  # "FILE", gsub(".ptx",".asc",input_file),
                  "-CLEAR",
                  sep = " "))
  
  while (term==1) sys.sleep(10)
  
  dir.create( substr(input_file, 0,nchar(input_file)-4))
  file.move(list.files("input", pattern = gsub("input/","",substr(input_file, 0,nchar(input_file)-4)), full.names = TRUE)[-1],
            substr(input_file, 0,nchar(input_file)-4))
  
}


readTLS<-function(output_file){
  
  dat<-data.table::fread(output_file, skip= 11, header = FALSE)
  colnames(dat)[1:7]<-c("X","Y","Z", "R","G","B","I")
  
  return(dat)
  
}

readTLSnorms<-function(output_file, cols_num){
  
  # dat<-data.table::fread(output_file, header = FALSE)

if(cols_num==7){
  dat<-vroom(output_file, delim = " ",
             col_names = c("X","Y","Z", "R","G","B","I","nX","nY","nZ"),
             # col_types= cols(),
             col_types= c(X='d',Y='d',Z='d', R='i',G='i',B='i',I='d',
                          nX='d',nY='d',nZ='d'),
             progress=FALSE)
} else {
  dat<-vroom(output_file, delim = " ",
             col_names = c("X","Y","Z","I","nX","nY","nZ"),
             # col_types= cols(),
             col_types= c(X='d',Y='d',Z='d',I='d',
                          nX='d',nY='d',nZ='d'),
             progress=FALSE)
}
# colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")

return(dat)

}

angleCalc<-function(dat, center, scatterLim=85, cols_num){
  #column naming
  # dat<-dat.out<-NULL
  # dat<-data.table::fread(output_file, header = FALSE)
  
  gc()
  # dat<-vroom(output_file, delim = " ", 
  #            col_names = c("X","Y","Z", "R","G","B","I","nX","nY","nZ"), 
  #            # col_types= cols(),
  #            col_types= c(X='d',Y='d',Z='d', R='i',G='i',B='i',I='d',
  #                         nX='d',nY='d',nZ='d'),
  #            progress=FALSE)
  # colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")
  
  # dat$dipDir_deg <- (pi+dat$dipDir_rad) * (180/pi)
  # dat$dip_deg <- (pi+dat$dip_rad) * (180/pi)
  
  # hist(dat$dipDir_deg)
  # hist(dat$dip_deg)
  
  #calculate the radius, zenith, and azimuth angles from XYZ coordinates
  dat$inc<-RZA(dat[,1:3])
  
  # estimate scattering angle and filter, removing steep angles
  dat$scatter<-scatter(dat, cols_num)
  dat<-dat[dat$scatter<=scatterLim,]
  dat<-na.omit(dat)
  
  #Convert normals to leaf orientation and leaf angle
  dat.out<-cbind(dat[,1:3],
                 dip=ConvertNormalToDipAndDipDir(dat[,c('nX','nY','nZ')])[,2])
  
  gc()
  
  dat.out<-dat.out[!is.nan(dat.out$dip),]
  return(dat.out)
  remove(dat, dat.out)
  
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
                  "-CLEAR",
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
                  "-CLEAR",
                  sep = " "))  
  
  while (term==1) sys.sleep(10)
  
}

rfPrep<-function(c2c.file){
  
  dat.temp<-dat.temp1<-dat.temp2<-dat.temp3<-NULL
  gc()
  dat.temp1<-data.table::fread(gsub(".asc","_0_10_NORM.asc",c2c.file),header = FALSE)
  # dat.temp1<-readr::read_table(gsub(".asc","_0_10_NORM.asc",c2c.file))
  # dat.temp1<-vroom(gsub(".asc","_0_10_NORM.asc",c2c.file),
  #                  col_names = c("X","Y","Z","angle","nX10","nY10","nZ10"),
  #                  col_types = c(X='d',Y='d',Z='d',angle='d',
  #                                nX10='d',nY10='d',nZ10='d'),
  #                  altrep = FALSE)
  colnames(dat.temp1)[1:7]<-c("X","Y","Z","angle","nX10","nY10","nZ10")
  
  # colnames(dat.temp1)[1:7]<-c(X='d',Y='d',Z='d',angle='d',
  #                             nX10='d',nY10='d',nZ10='d')
  
  dat.temp2<-data.table::fread(gsub(".asc","_0_50_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
  # dat.temp2<-vroom(gsub(".asc","_0_50_NORM.asc",c2c.file), 
  #                  col_select = c(5,6,7),
  #                  col_names = c("X","Y","Z","angle","nX50","nY50","nZ50"),
  #                  col_types = c(nX50='d',nY50='d',nZ50='d'),
  #                  altrep = FALSE)
  colnames(dat.temp2)[1:3]<-c("nX50","nY50","nZ50")
  
  dat.temp3<-data.table::fread(gsub(".asc","_0_75_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
  # dat.temp3<-vroom(gsub(".asc","_0_75_NORM.asc",c2c.file), 
  #                  col_select = c(5,6,7),
  #                  col_names = c("X","Y","Z","angle","nX75","nY75","nZ75"),
  #                  col_types = c(nX75='d',nY75='d',nZ75='d'),
  #                  altrep = FALSE)
  colnames(dat.temp3)[1:3]<-c("nX75","nY75","nZ75")
  
  
  dat.temp<-as.data.frame(cbind(dat.temp1,dat.temp2,dat.temp3))
  remove(dat.temp1,dat.temp2,dat.temp3)
  gc()
  
  # return(dat.temp)
  data.table::fwrite(dat.temp, file = gsub(".asc","_rf_prep.asc",c2c.file), sep = " ", row.names = FALSE)
  remove(dat.temp)
  gc()
  
  file.remove(gsub(".asc","_0_10_NORM.asc",c2c.file),
              gsub(".asc","_0_50_NORM.asc",c2c.file),
              gsub(".asc","_0_75_NORM.asc",c2c.file))
  
}


angleCalcWrite<-function(output_file, center, scatterLim, angle.file.name){
  dat.angle<-NULL
  dat.angle<-angleCalc(output_file,
                       center,
                       scatterLim)
  
  fwrite(dat.angle, file = angle.file.name, 
         sep = " ", row.names = FALSE)
  
  remove(dat.angle)
}

rf_predict<-function(c2c.file, 
                     rf_model,
                     rf_model_path,
                     class.file.name){
  dat.rf<-dat.rf.out<-NULL
  gc()
  
  if(is.null(rf_model)){
    rf_model<-NULL
    rf_model<-readRDS(rf_model_path)
  }
  
  
  # dat.rf<-as.data.frame(na.omit(rfPrep(c2c.file)))
  dat.rf<-as.data.frame(na.omit(fread(gsub(".asc","_rf_prep.asc",c2c.file))))
  dat.rf$predict<-predict(rf_model, dat.rf)
  
  # rf_model<-NULL
  gc()
  
  dat.rf.out<-cbind(dat.rf[,1:4], dat.rf$predict)
  colnames(dat.rf.out)[4:5]<-c("dip_deg", "class")
  
  remove(dat.rf)
  gc()
  
  data.table::fwrite(dat.rf.out, file = class.file.name, sep = " ", row.names = FALSE)
  remove(dat.rf.out)
  dat.rf<-dat.rf.out<-NULL
  
  # file.remove(gsub(".asc","_rf_prep.asc",c2c.file))
  
  gc()
}

voxel_beta_fit<-function(class.file.name, 
                         voxRes=0.1, 
                         minVoxDensity=5, 
                         correct.topography=TRUE){
  dat.class.in<-voxels<-LAD<-LAD_lim<-m<-beta<-NULL
  gc()
  dat.class.in<-vroom::vroom(class.file.name, delim = " ", 
                             col_names = c("X","Y","Z", "dip_deg","class"), 
                             col_types= c(X='d',Y='d',Z='d',dip_deg='d',class='i'),
                             progress=FALSE, skip=1,
                             altrep = FALSE)
  
  # if(file.exists(class.file.name)) dat.class<-vroom(class.file.name)
  dat.class<-NULL
  print("Correct topography...")
  #Correct topography and calculate the LAD and vertical LAD
  if(correct.topography) dat.class<-normalize_topography(dat.class.in) else dat.class<-dat.class.in
  print("Done")
  # fwrite(dat.class,  file = gsub(".asc", "_topo_correct.asc",class.file.name),
  #        sep = " ", row.names = FALSE)
  
  
  print("Voxelize angle estimates...")
  # leaf angle voxelation and density normalization
  
  if(correct.topography) dat.class.ag<- dat.class[dat.class$Z>1,] else dat.class.ag<-dat.class
  
  voxels<-LAvoxel(na.omit(dat.class.ag), voxRes)
  
  voxels<-voxels[voxels$n>quantile(voxels$n,0.01),]
  voxels<-voxels[voxels$n>minVoxDensity,]
  
  print("Done")
  
  fwrite(voxels,  file = gsub(".asc", "_voxels.asc",class.file.name),
         sep = " ", row.names = FALSE)
  
  #simulate LAD from voxel statistics
  print("Simulate LAD from voxel statistics...")
  
  if(nrow(voxels)>0){
    LAD<-sim_LAD(voxels)
    LAD<-na.exclude(LAD[LAD$a>=0 & LAD$a<=90,])
    fwrite(LAD,  file = gsub(".asc", "_LAD.txt",class.file.name),
           sep = " ", row.names = FALSE)
    
    LAD_lim<-data.frame(a=LAD$a)
    
    #fit beta function and get beta parameters from LAD
    #FIT BETA DISTRIBUTION
    
    # get_beta<-function(x){
    m<-NULL
    print("Fit Beta distribution to LAD...")
    # if (nrow(LAD_lim)>1){
    m <- fitdistrplus::fitdist(as.numeric(LAD_lim$a)/90, 'beta', method='mle')
    
    # Get alpha and beta parametrs
    alpha0 <- m$estimate[1] # parameter 1
    beta0 <- m$estimate[2] # parameter 2
    beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01),
                                                     alpha0,beta0))
    param<-data.frame(alpha0,beta0)
    print("Done")
    
    fwrite(beta,  file = gsub(".asc", paste("_Beta_distribution_alpha_", 
                                            round(param[1],2),"_beta_",round(param[2],2),
                                            ".txt", sep=""),class.file.name),
           sep = " ", row.names = FALSE)
    
    remove(m)
    gc()
    # } 
  }
  
  # gc()
}



TLSLeAF<-function(input_file,
                  overwrite=TRUE,
                  center, 
                  scatterLim=85,
                  SS=0.02, 
                  scales=c(0.1,0.5,0.75),
                  rf_model=NULL,
                  rf_model_path=NULL,
                  correct.topography=TRUE,
                  voxRes = 0.1,
                  minVoxDensity=5,
                  superDF=TRUE,
                  clean=TRUE,...){
  
  # keep.vars<-ls()
  
  setDTthreads(threads = 1)
  # getDTthreads(verbose=TRUE)
  
  #names of output files
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  gc()
  
  cols_num<-length(fread(input_file, skip=11, nrows=1, header=FALSE))
  
  print("Calculate normals for gridded TLS point cloud...")
  #Calculate normals for gridded TLS point cloud
  if(!file.exists(output_file)|
     overwrite) normalCalc(input_file)
  print("Done")
  
  while(!file.exists(output_file)) Sys.sleep(10)
  while(file.size(output_file)<file.size(input_file)) Sys.sleep(10)
  while((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
    rstudioapi::terminalBusy(rstudioapi::terminalList())])>0)) Sys.sleep(1)
  
  if((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
    rstudioapi::terminalBusy(rstudioapi::terminalList())])==0)) rstudioapi::terminalKill(rstudioapi::terminalBusy(rstudioapi::terminalList()))
  
  print("Calculate scattering angle and leaf angle...")
  #Calculate scattering angle and leaf angle
  if((!file.exists(angle.file.name)|
     overwrite)){
    # remove(dat)
    # gc()
    
    dat<-readTLSnorms(output_file, cols_num)
    
    dat.angle<-angleCalc(dat,
                         center,
                         scatterLim,
                         cols_num)
    
    fwrite(dat.angle, file = angle.file.name, 
           sep = " ", row.names = FALSE)
    
    remove(dat,dat.angle)
    gc()
    
    # print("Calculating surface angles...")
    # angleCalcWrite(output_file, center, scatterLim, angle.file.name)
    
    
    
  }
  print("Done")
  
  while(!file.exists(angle.file.name)) Sys.sleep(10)
  
  print("Calculating normals for leaf and wood classification...")
  #Classify wood and leaf from random forest classifier
  if(file.exists(angle.file.name)&
     !(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))|
     overwrite) classMetricCalc(angle.file.name, SS, scales)
  
  print("Done")
  
  while(!file.exists(c2c.file)) Sys.sleep(10)
  while(!(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
          file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
          file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))) Sys.sleep(10)
  
  if((file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))&
     !file.exists(gsub(".asc","_rf_prep.asc",c2c.file))|
     overwrite) {
    
    print("Prepping for leaf and wood classification...")
    rfPrep(c2c.file)
    
    print("Done")  
  }
  
  while(!file.exists(gsub(".asc","_rf_prep.asc",c2c.file))) Sys.sleep(10)
  
  if((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
    rstudioapi::terminalBusy(rstudioapi::terminalList())])==0)) rstudioapi::terminalKill(rstudioapi::terminalBusy(rstudioapi::terminalList()))
  
  if(!file.exists(class.file.name)&
     file.exists(gsub(".asc","_rf_prep.asc",c2c.file))|
     overwrite) {
    
    print("Classifying leaf and wood points...")
    rf_predict(c2c.file=c2c.file,
               rf_model,
               rf_model_path = rf_model_path,
               class.file.name=class.file.name)
    
    
  }
  
  status<-TRUE
  # else 
  print("Done")
  
  while(!file.exists(class.file.name)) Sys.sleep(5)
  
  print("Voxelizing and calculating LAD...")
  if(file.exists(class.file.name)&status){
    
    voxel_beta_fit(class.file.name, 
                   voxRes, 
                   minVoxDensity, 
                   correct.topography)
  }
  print("Done")
  
  
  print("Calculating G-function...")
  if(file.exists(class.file.name)&
     file.exists(gsub(".asc", "_LAD.txt",class.file.name))&
     file.exists(voxels=gsub(".asc", "_voxels.asc",class.file.name))){
    G_calculations(dat=class.file.name,
                   LAD=gsub(".asc", "_LAD.txt",class.file.name),
                   voxels=gsub(".asc", "_voxels.asc",class.file.name))
  }
  
  # while(!file.exists(gsub(".asc", paste("_Beta_distribution_alpha_",
  #                                       round(param[1],2),"_beta_",round(param[2],2),
  #                                       ".txt", sep=""),class.file.name))) Sys.sleep(5)
  
  if (superDF){
    
    dat.class<-LAD<-voxels<-G<-m<-beta<-NULL
    
    dat.class<-fread(class.file.name)
    LAD<-fread(gsub(".asc", "_LAD.txt",class.file.name))
    voxels<-fread(gsub(".asc", "_voxels.asc",class.file.name))
    G<-fread(gsub("angles_class_voxels.asc", "G_function.txt",gsub(".asc", "_voxels.asc",class.file.name)))
    
    # get_beta<-function(x){
    m<-NULL
    print("Fit Beta distribution to LAD...")
    # if (nrow(LAD_lim)>1){
    m <- fitdistrplus::fitdist(as.numeric(LAD$a)/90, 'beta', method='mle')
    
    # Get alpha and beta parametrs
    alpha0 <- m$estimate[1] # parameter 1
    beta0 <- m$estimate[2] # parameter 2
    beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01),
                                                     alpha0,beta0))
    param<-data.frame(alpha0,beta0)
    
    print("Finalizing TLSLeAF object...")
    
    TLSLeAF.dat<-new("TLSLeAF",
                     parameters=data.frame(c(file=input_file,
                                             center,
                                             scatterLim=85,
                                             SS=0.02,
                                             scales=c(0.1,0.5,0.75),
                                             voxRes=voxRes,
                                             superDF=TRUE)),
                     dat=as.data.frame(dat.class),
                     voxels=as.data.frame(voxels),
                     LAD=LAD,
                     Beta_parameters=param,
                     beta=beta,
                     G=G)
    
    # remove(dat)
    remove(dat.class)
    remove(voxels)
    remove(LAD)
    gc()
    return(TLSLeAF.dat)
  }
  
  if(clean){
    clean.temp(output_file, c2c.file, clean=TRUE)
  }
  
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

clean.temp<-function(output_file, c2c.file, clean=TRUE){
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


pgapF<-function (tls.scan.sub, row_res, col_res, z_res) {
  total<-list<-pgap<-NULL
  total<-(row_res/150)*col_res*5
  list <- seq(floor(min(tls.scan.sub[3])), 
              ceiling(max(tls.scan.sub[3])), 
              by = z_res)
  
  pgap <- NULL
  for (i in 1:length(list)) {
    pgap[i] <- 1-(length(subset(tls.scan.sub$Z, tls.scan.sub$Z<list[i])))/total
  }
  pgap<- (1-max(pgap))+pgap
  return(pgap)
}

a_vai<-function(tls.scan, a_ls, z_res, col_res, row_res){
  vai_a<-list()
  for (a in 1:length(a_ls)){
    
    tls.scan.sub<-subset(tls.scan,tls.scan$inc>a_ls[a]&tls.scan$inc<(a_ls[a]+5))
    list <- seq(floor(min(tls.scan.sub[3])), 
                ceiling(max(tls.scan.sub[3])), 
                by = z_res)
    
    pgap<-pgapF(tls.scan.sub, row_res, col_res, z_res)
    
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
      # tls.scan = gsub(".ptx","", files[l]),
      a = a_ls[a],
      z = f.prime.sub[,2],
      dx = f.prime.sub[,1])
    
    
  }
  return(data.frame(do.call(rbind,vai_a)))
}

vai<-function (tls.scan, a_range, a_bin, z_res, col_res, row_res, name){
  # tls.scan<-RZA(tls.scan)
  a_ls<-seq(min(a_range),max(a_range), by = a_bin) # zenith angle range
  
  vai_out<-data.frame(a_vai(tls.scan, a_ls, z_res, col_res, row_res),
                      tls.scan = name)
  return(vai_out)
}



vai_weighted <- function(a_vai, adj_z, z_res) {
  a_vai <- a_vai[a_vai$z>0,]
  a_vai$z_adj<-a_vai$z+adj_z
  a_vai$z_bin<-floor(a_vai$z_adj)
  #weight the profiles by zenith angle
  a_vai$weight<-sin(a_vai$a*(pi/180))
  a_vai_weight<-a_vai %>%
    group_by(tls.scan, z_bin) %>% 
    mutate(dxw = weighted.mean(dx, weight))
  return(a_vai_weight)
}


ptx.header<-function(input_file){
  list(
    name = input_file,
    col_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, header = FALSE)[1,]),
    row_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, header = FALSE)[2,]),
    scan_center = as.matrix(read.csv(input_file, sep = "", skip = 2,nrows = 1, header = FALSE)),
    reg = as.matrix(read.csv(input_file, sep = "", skip = 3,nrows = 3, header = FALSE)),
    trans_matrix = as.matrix(read.csv(input_file, sep = "", skip = 6,nrows = 4, header = FALSE))
  )
}

pgap.angle<-function(tls.scan,
                     z_res=1,
                     adj_z=1.6,
                     a_range=c(10,60),
                     a_bin = 5,
                     plane = FALSE,
                     header,...){
  
  
  
  #calculate inclination angle
  inc<-RZA(tls.scan, deg=TRUE)
  tls.scan$inc<-inc
  remove(inc)
  gc()
  gc()
  
  # zenith angle range
  a_ls<-seq(min(a_range),max(a_range), by = a_bin) 
  
  pgap_ls<-list()
  
  pgap_ls<-lapply(a_ls, function(x){
    tls.scan.sub<-subset(tls.scan,tls.scan$inc>x&tls.scan$inc<(x+a_bin))
    tls.scan.sub<-na.exclude(tls.scan.sub)
    # if(nrow(tls.scan.sub) ==0) next
    list <- seq(floor(min(tls.scan.sub[3])), 
                ceiling(max(tls.scan.sub[3])), 
                by = z_res)
    
    pgap<-pgapF(tls.scan.sub, header$row_res, header$col_res, z_res)
    # plot(pgap)
    
    remove(tls.scan.sub)
    gc()
    
    return(data.frame(list,pgap, inc_bin=x))
  })
  remove(tls.scan)
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


G_calculations<-function(dat,
                         LAD,
                         voxels){
  
  G_ls<-list()
  if(file.exists(dat)){
    # dat<-NULL
    gc()
    dat<-fread(dat)
    dat<-dat[dat$class==0,]
    
    if(nrow(dat)>0){
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
      for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(dip_bin=G_test$dip_bin[G_test$dip_bin==c(0:90)[i]],
                                                         count=G_test$count[G_test$dip_bin==c(0:90)[i]],
                                                         p=G_test$p[G_test$dip_bin==c(0:90)[i]],
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
      
      # plot(G_calc$inc_bin,
      #      G_calc$G_p)
      
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
    } 
  }
  
  if(file.exists(voxels)){
    dat<-fread(voxels)
    if(nrow(dat)>0){
      
      dat$inc<-RZA(dat[,1:3], deg = TRUE)
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
      for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(dip_bin=G_test$dip_bin[G_test$dip_bin==c(0:90)[i]],
                                                         count=G_test$count[G_test$dip_bin==c(0:90)[i]],
                                                         p=G_test$p[G_test$dip_bin==c(0:90)[i]],
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
    }
  }
  
  fwrite(do.call(rbind, G_ls),
         file = gsub("angles_class_voxels.asc", "G_function.txt",voxels),
         sep= " ")
  dat<-NULL
  G_ls<-NULL
  gc()
  # return(do.call(rbind, G_ls))
}
=======
ConvertNormalToDipAndDipDir<-function(x) {
  
  Nsign<-dipDir_rad<-dip_rad<-dip_deg<-NULL
  
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
  r<-inc<-NULL
  
  colnames(dat)[1:3] <- c("X","Y","Z")
  r <- sqrt(dat$X^2 + dat$Y^2 +dat$Z^2)
  
  inc<-matrix(NA, nrow=length(r))
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

scatter<-function(dat, cols_num){
  my.dt<-NULL
  # time<-Sys.time()
  # dat<-cbind(dat, )
  dat.add<-dat[,1:3]/sqrt(rowSums(dat[,1:3]^2))
  colnames(dat.add)<-c("sX","sY","sZ")
  dat<-cbind(dat,dat.add)
  remove(dat.add)
  
  dat$dot<-geometry::dot(dat[,cols_num+1:3],dat[,c('nX','nY','nZ')],d=2)
  my.dt <- as.data.table(cbind(abs(dat$dot),1))
  dat$dot<-my.dt[,row.min:=pmin(V1,V2)]$row.min
  return(acos(dat$dot) * (180/pi))
  remove(my.dt)
  # print(Sys.time()-time)#
}

normalize_topography<-function(x, res = 5){
  las<-r<-topo<-topo.df<-ground<-topo.las.r<-r.p<-topo.p<-slope<-poly.id<-NULL
  
  # defaultW <- getOption("warn")
  # options(warn = -1)
  # crs_tls <- sp::CRS("+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs")
  crs_tls<-sp::CRS("+init=epsg:26918")
  # crs_tls
  colnames(x)[1:3]<-c("X","Y","Z")
  las<-LAS(x[,1:3], proj4string = crs_tls)
  # projection(las)<-crs
  las@data$dip_deg<-x[,4]
  las@data$class<-x[,5]
  
  r <- raster(xmn=-200, xmx=200, ymn=-200, ymx=200, resolution = res)
  topo<-grid_metrics(las, quantile(Z, 0.01), r)
  crs(topo)<-crs_tls
  # plot(topo)
  
  topo.df<-as.data.frame(rasterToPoints(topo))
  colnames(topo.df)<-c("X","Y","Z")
  
  ws <- seq(3,12, 3)
  th <- seq(0.1, 1.5, length.out = length(ws))
  
  topo.las<-LAS(topo.df)
  crs(topo.las)<-crs_tls
  
  ground<-classify_ground(topo.las, pmf(ws, th), last_returns = FALSE)
  topo.las.r<-grid_terrain(ground, r, kriging(k = 10))
  
  # topo<-focal(topo, matrix(1,3,3))
  # topo.smooth<-disaggregate(topo, 5, method='bilinear')
  # # density.r<-grid_density(las, topo.smooth)
  # # plot(topo, col = viridis::viridis(250))
  # plot(topo.las.r)
  
  r.p<-buffer(topo.las.r, width=1)
  r.p[!is.na(r.p)]<-1
  
  topo.las.p<-rasterToPolygons(r.p, dissolve = TRUE)
  # lines(topo.las.p)
  
  las<-clip_polygon(las, topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,1], topo.las.p@polygons[[1]]@Polygons[[1]]@coords[,2])
  
  # plot(las)
  las<- las - topo.las.r
  # plot(las)
  
  slope<-terrain(topo, opt = "slope", unit = "degrees", neighbors = 8)
  # plot(slope)
  
  topo[is.na(slope)]<-NA
  topo.p<-topo
  topo.p[!is.na(topo.p)]<-1
  
  topo.las.p<-rasterToPolygons(topo.p, dissolve = TRUE)
  # lines(topo.las.p)
  
  poly.id<-length(topo.las.p@polygons[[1]]@Polygons)
  
  topo.las.p@polygons[[1]]@Polygons
  
  a.ls<-c(do.call(rbind, lapply(1:poly.id, function(x){
    area.sub<-topo.las.p@polygons[[1]]@Polygons[[x]]
    return(area.sub@area)
  })))
  
  max.poly.id<-(1:poly.id)[max(a.ls)==a.ls]
  
  las<-clip_polygon(las, topo.las.p@polygons[[1]]@Polygons[[max.poly.id]]@coords[,1], topo.las.p@polygons[[1]]@Polygons[[max.poly.id]]@coords[,2])
  # plot(las)
  return(las@data)
  
  # remove(las)
  
}

LAvoxel<-function(x,res = 0.1){
  # defaultW <- getOption("warn")
  # options(warn = -1)
  las<-voxels<-NULL
  colnames(x)[1:3]<-c("X","Y","Z")
  x<-x[x$class==0,]
  crs_tls<-sp::CRS("+init=epsg:26918")
  las<-LAS(x[,1:3],proj4string = crs_tls)
  
  # crs(las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  # las@data$Z<-x$z_cor
  las@data$dip_deg<-x$dip_deg
  las@data$count<-1
  
  voxels<-voxel_metrics(las, ~mean(dip_deg), res=res)
  colnames(voxels)[4]<-"dip_dir"
  voxels$dip_dir_sd<-voxel_metrics(las, ~sd(dip_deg), res=res)[,4]
  voxels$n<-voxel_metrics(las, ~sum(count), res=res)[,4]
  
  # options(warn = defaultW)
  return(voxels)
  
  remove(las)
  remove(voxels)
}

sim_LAD<-function(x){
  voxels.dup<-NULL
  vox_ls<-list()
  for(i in 1:10){
    voxels.dup<-x
    voxels.dup$a<-apply(voxels.dup,1,FUN = function(xx)  rnorm(1, mean=xx[4], sd=xx[5]))
    vox_ls[[i]]<-voxels.dup
    remove(voxels.dup)
  }
  return(do.call(rbind, vox_ls))
  remove(vox_ls)
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
                  "-AUTO_SAVE OFF",
                  "-COMPUTE_NORMALS",
                  "-O", input_file, #open the subsampled file
                  # "-NORMALS_TO_DIP",
                  "-SAVE_CLOUDS", "FILE", gsub(".ptx",".asc",input_file),
                  "-CLEAR",
                  sep = " "))
  
  while (term==1) sys.sleep(10)
}


readTLS<-function(output_file){
  
  dat<-data.table::fread(output_file, skip= 11, header = FALSE)
  colnames(dat)[1:7]<-c("X","Y","Z", "R","G","B","I")
  
  return(dat)
  
}

readTLSnorms<-function(output_file, cols_num){
  
  # dat<-data.table::fread(output_file, header = FALSE)

if(cols_num==7){
  dat<-vroom(output_file, delim = " ",
             col_names = c("X","Y","Z", "R","G","B","I","nX","nY","nZ"),
             # col_types= cols(),
             col_types= c(X='d',Y='d',Z='d', R='i',G='i',B='i',I='d',
                          nX='d',nY='d',nZ='d'),
             progress=FALSE)
} else {
  dat<-vroom(output_file, delim = " ",
             col_names = c("X","Y","Z","I","nX","nY","nZ"),
             # col_types= cols(),
             col_types= c(X='d',Y='d',Z='d',I='d',
                          nX='d',nY='d',nZ='d'),
             progress=FALSE)
}
# colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")

return(dat)

}

angleCalc<-function(dat, center, scatterLim=85, cols_num){
  #column naming
  # dat<-dat.out<-NULL
  # dat<-data.table::fread(output_file, header = FALSE)
  
  gc()
  # dat<-vroom(output_file, delim = " ", 
  #            col_names = c("X","Y","Z", "R","G","B","I","nX","nY","nZ"), 
  #            # col_types= cols(),
  #            col_types= c(X='d',Y='d',Z='d', R='i',G='i',B='i',I='d',
  #                         nX='d',nY='d',nZ='d'),
  #            progress=FALSE)
  # colnames(dat)[1:10]<-c("X","Y","Z", "R","G","B","I","nX","nY","nZ")
  
  # dat$dipDir_deg <- (pi+dat$dipDir_rad) * (180/pi)
  # dat$dip_deg <- (pi+dat$dip_rad) * (180/pi)
  
  # hist(dat$dipDir_deg)
  # hist(dat$dip_deg)
  
  #calculate the radius, zenith, and azimuth angles from XYZ coordinates
  dat$inc<-RZA(dat[,1:3])
  
  # estimate scattering angle and filter, removing steep angles
  dat$scatter<-scatter(dat, cols_num)
  dat<-dat[dat$scatter<=scatterLim,]
  dat<-na.omit(dat)
  
  #Convert normals to leaf orientation and leaf angle
  dat.out<-cbind(dat[,1:3],
                 dip=ConvertNormalToDipAndDipDir(dat[,c('nX','nY','nZ')])[,2])
  
  gc()
  
  dat.out<-dat.out[!is.nan(dat.out$dip),]
  return(dat.out)
  remove(dat, dat.out)
  
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
                  "-CLEAR",
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
                  "-CLEAR",
                  sep = " "))  
  
  while (term==1) sys.sleep(10)
  
}

rfPrep<-function(c2c.file){
  
  dat.temp<-dat.temp1<-dat.temp2<-dat.temp3<-NULL
  gc()
  dat.temp1<-data.table::fread(gsub(".asc","_0_10_NORM.asc",c2c.file),header = FALSE)
  # dat.temp1<-readr::read_table(gsub(".asc","_0_10_NORM.asc",c2c.file))
  # dat.temp1<-vroom(gsub(".asc","_0_10_NORM.asc",c2c.file),
  #                  col_names = c("X","Y","Z","angle","nX10","nY10","nZ10"),
  #                  col_types = c(X='d',Y='d',Z='d',angle='d',
  #                                nX10='d',nY10='d',nZ10='d'),
  #                  altrep = FALSE)
  colnames(dat.temp1)[1:7]<-c("X","Y","Z","angle","nX10","nY10","nZ10")
  
  # colnames(dat.temp1)[1:7]<-c(X='d',Y='d',Z='d',angle='d',
  #                             nX10='d',nY10='d',nZ10='d')
  
  dat.temp2<-data.table::fread(gsub(".asc","_0_50_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
  # dat.temp2<-vroom(gsub(".asc","_0_50_NORM.asc",c2c.file), 
  #                  col_select = c(5,6,7),
  #                  col_names = c("X","Y","Z","angle","nX50","nY50","nZ50"),
  #                  col_types = c(nX50='d',nY50='d',nZ50='d'),
  #                  altrep = FALSE)
  colnames(dat.temp2)[1:3]<-c("nX50","nY50","nZ50")
  
  dat.temp3<-data.table::fread(gsub(".asc","_0_75_NORM.asc",c2c.file), select = c(5,6,7),header = FALSE)
  # dat.temp3<-vroom(gsub(".asc","_0_75_NORM.asc",c2c.file), 
  #                  col_select = c(5,6,7),
  #                  col_names = c("X","Y","Z","angle","nX75","nY75","nZ75"),
  #                  col_types = c(nX75='d',nY75='d',nZ75='d'),
  #                  altrep = FALSE)
  colnames(dat.temp3)[1:3]<-c("nX75","nY75","nZ75")
  
  
  dat.temp<-as.data.frame(cbind(dat.temp1,dat.temp2,dat.temp3))
  remove(dat.temp1,dat.temp2,dat.temp3)
  gc()
  
  # return(dat.temp)
  data.table::fwrite(dat.temp, file = gsub(".asc","_rf_prep.asc",c2c.file), sep = " ", row.names = FALSE)
  remove(dat.temp)
  gc()
  
  file.remove(gsub(".asc","_0_10_NORM.asc",c2c.file),
              gsub(".asc","_0_50_NORM.asc",c2c.file),
              gsub(".asc","_0_75_NORM.asc",c2c.file))
  
}


angleCalcWrite<-function(output_file, center, scatterLim, angle.file.name){
  dat.angle<-NULL
  dat.angle<-angleCalc(output_file,
                       center,
                       scatterLim)
  
  fwrite(dat.angle, file = angle.file.name, 
         sep = " ", row.names = FALSE)
  
  remove(dat.angle)
}

rf_predict<-function(c2c.file, 
                     rf_model,
                     rf_model_path,
                     class.file.name){
  dat.rf<-dat.rf.out<-NULL
  gc()
  
  if(is.null(rf_model)){
    rf_model<-NULL
    rf_model<-readRDS(rf_model_path)
  }
  
  
  # dat.rf<-as.data.frame(na.omit(rfPrep(c2c.file)))
  dat.rf<-as.data.frame(na.omit(fread(gsub(".asc","_rf_prep.asc",c2c.file))))
  dat.rf$predict<-predict(rf_model, dat.rf)
  
  # rf_model<-NULL
  gc()
  
  dat.rf.out<-cbind(dat.rf[,1:4], dat.rf$predict)
  colnames(dat.rf.out)[4:5]<-c("dip_deg", "class")
  
  remove(dat.rf)
  gc()
  
  data.table::fwrite(dat.rf.out, file = class.file.name, sep = " ", row.names = FALSE)
  remove(dat.rf.out)
  dat.rf<-dat.rf.out<-NULL
  
  # file.remove(gsub(".asc","_rf_prep.asc",c2c.file))
  
  gc()
}

voxel_beta_fit<-function(class.file.name, 
                         voxRes=0.1, 
                         minVoxDensity=5, 
                         correct.topography=TRUE){
  dat.class.in<-voxels<-LAD<-LAD_lim<-m<-beta<-NULL
  gc()
  dat.class.in<-vroom::vroom(class.file.name, delim = " ", 
                             col_names = c("X","Y","Z", "dip_deg","class"), 
                             col_types= c(X='d',Y='d',Z='d',dip_deg='d',class='i'),
                             progress=FALSE, skip=1,
                             altrep = FALSE)
  
  # if(file.exists(class.file.name)) dat.class<-vroom(class.file.name)
  dat.class<-NULL
  print("Correct topography...")
  #Correct topography and calculate the LAD and vertical LAD
  if(correct.topography) dat.class<-normalize_topography(dat.class.in) else dat.class<-dat.class.in
  print("Done")
  # fwrite(dat.class,  file = gsub(".asc", "_topo_correct.asc",class.file.name),
  #        sep = " ", row.names = FALSE)
  
  
  print("Voxelize angle estimates...")
  # leaf angle voxelation and density normalization
  
  if(correct.topography) dat.class.ag<- dat.class[dat.class$Z>1,] else dat.class.ag<-dat.class
  
  voxels<-LAvoxel(na.omit(dat.class.ag), voxRes)
  
  voxels<-voxels[voxels$n>quantile(voxels$n,0.01),]
  voxels<-voxels[voxels$n>minVoxDensity,]
  
  print("Done")
  
  fwrite(voxels,  file = gsub(".asc", "_voxels.asc",class.file.name),
         sep = " ", row.names = FALSE)
  
  #simulate LAD from voxel statistics
  print("Simulate LAD from voxel statistics...")
  
  if(nrow(voxels)>0){
    LAD<-sim_LAD(voxels)
    LAD<-na.exclude(LAD[LAD$a>=0 & LAD$a<=90,])
    fwrite(LAD,  file = gsub(".asc", "_LAD.txt",class.file.name),
           sep = " ", row.names = FALSE)
    
    LAD_lim<-data.frame(a=LAD$a)
    
    #fit beta function and get beta parameters from LAD
    #FIT BETA DISTRIBUTION
    
    # get_beta<-function(x){
    m<-NULL
    print("Fit Beta distribution to LAD...")
    # if (nrow(LAD_lim)>1){
    m <- fitdistrplus::fitdist(as.numeric(LAD_lim$a)/90, 'beta', method='mle')
    
    # Get alpha and beta parametrs
    alpha0 <- m$estimate[1] # parameter 1
    beta0 <- m$estimate[2] # parameter 2
    beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01),
                                                     alpha0,beta0))
    param<-data.frame(alpha0,beta0)
    print("Done")
    
    fwrite(beta,  file = gsub(".asc", paste("_Beta_distribution_alpha_", 
                                            round(param[1],2),"_beta_",round(param[2],2),
                                            ".txt", sep=""),class.file.name),
           sep = " ", row.names = FALSE)
    
    remove(m)
    gc()
    # } 
  }
  
  # gc()
}



TLSLeAF<-function(input_file,
                  overwrite=TRUE,
                  center, 
                  scatterLim=85,
                  SS=0.02, 
                  scales=c(0.1,0.5,0.75),
                  rf_model=NULL,
                  rf_model_path=NULL,
                  correct.topography=TRUE,
                  voxRes = 0.1,
                  minVoxDensity=5,
                  superDF=TRUE,
                  clean=TRUE,
                  OS='mac',...){
  
  # keep.vars<-ls()
  
  setDTthreads(threads = 1)
  # getDTthreads(verbose=TRUE)
  
  #names of output files
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  gc()
  
  cols_num<-length(fread(input_file, skip=11, nrows=1, header=FALSE))
  
  print("Calculate normals for gridded TLS point cloud...")
  #Calculate normals for gridded TLS point cloud
  if(!file.exists(output_file)|
     overwrite) normalCalc(input_file)
  print("Done")
  
  while(!file.exists(output_file)) Sys.sleep(10)
  while(file.size(output_file)<file.size(input_file)) Sys.sleep(10)
  
  if(OS=='mac'){
    while((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
      rstudioapi::terminalBusy(rstudioapi::terminalList())])>0)) Sys.sleep(1)
    
    if((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
      rstudioapi::terminalBusy(rstudioapi::terminalList())])==0)) rstudioapi::terminalKill(rstudioapi::terminalBusy(rstudioapi::terminalList()))
    
  }
  
  print("Calculate scattering angle and leaf angle...")
  #Calculate scattering angle and leaf angle
  if((!file.exists(angle.file.name)|
     overwrite)){
    # remove(dat)
    # gc()
    
    dat<-readTLSnorms(output_file, cols_num)
    
    dat.angle<-angleCalc(dat,
                         center,
                         scatterLim,
                         cols_num)
    
    fwrite(dat.angle, file = angle.file.name, 
           sep = " ", row.names = FALSE)
    
    remove(dat,dat.angle)
    gc()
    
    # print("Calculating surface angles...")
    # angleCalcWrite(output_file, center, scatterLim, angle.file.name)
    
    
    
  }
  print("Done")
  
  while(!file.exists(angle.file.name)) Sys.sleep(10)
  
  print("Calculating normals for leaf and wood classification...")
  #Classify wood and leaf from random forest classifier
  if(file.exists(angle.file.name)&
     !(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
       file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))|
     overwrite) classMetricCalc(angle.file.name, SS, scales)
  
  print("Done")
  
  while(!file.exists(c2c.file)) Sys.sleep(10)
  while(!(file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
          file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
          file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))) Sys.sleep(10)
  
  if((file.exists(gsub(".asc","_0_10_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_50_NORM.asc",c2c.file))&
      file.exists(gsub(".asc","_0_75_NORM.asc",c2c.file)))&
     !file.exists(gsub(".asc","_rf_prep.asc",c2c.file))|
     overwrite) {
    
    print("Prepping for leaf and wood classification...")
    rfPrep(c2c.file)
    
    print("Done")  
  }
  
  while(!file.exists(gsub(".asc","_rf_prep.asc",c2c.file))) Sys.sleep(10)
  
  if(OS=='mac'){
    if((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
      rstudioapi::terminalBusy(rstudioapi::terminalList())])==0)) rstudioapi::terminalKill(rstudioapi::terminalBusy(rstudioapi::terminalList()))
  }
  
  if(!file.exists(class.file.name)&
     file.exists(gsub(".asc","_rf_prep.asc",c2c.file))|
     overwrite) {
    
    print("Classifying leaf and wood points...")
    rf_predict(c2c.file=c2c.file,
               rf_model,
               rf_model_path = rf_model_path,
               class.file.name=class.file.name)
    
    
  }
  
  status<-TRUE
  # else 
  print("Done")
  
  while(!file.exists(class.file.name)) Sys.sleep(5)
  
  print("Voxelizing and calculating LAD...")
  if(file.exists(class.file.name)&status){
    
    voxel_beta_fit(class.file.name, 
                   voxRes, 
                   minVoxDensity, 
                   correct.topography)
  }
  print("Done")
  
  
  print("Calculating G-function...")
  if(file.exists(class.file.name)&
     file.exists(gsub(".asc", "_LAD.txt",class.file.name))&
     file.exists(voxels=gsub(".asc", "_voxels.asc",class.file.name))){
    G_calculations(dat=class.file.name,
                   LAD=gsub(".asc", "_LAD.txt",class.file.name),
                   voxels=gsub(".asc", "_voxels.asc",class.file.name))
  }
  
  # while(!file.exists(gsub(".asc", paste("_Beta_distribution_alpha_",
  #                                       round(param[1],2),"_beta_",round(param[2],2),
  #                                       ".txt", sep=""),class.file.name))) Sys.sleep(5)
  
  if (superDF){
    
    dat.class<-LAD<-voxels<-G<-m<-beta<-NULL
    
    dat.class<-fread(class.file.name)
    LAD<-fread(gsub(".asc", "_LAD.txt",class.file.name))
    voxels<-fread(gsub(".asc", "_voxels.asc",class.file.name))
    G<-fread(gsub("angles_class_voxels.asc", "G_function.txt",gsub(".asc", "_voxels.asc",class.file.name)))
    
    # get_beta<-function(x){
    m<-NULL
    print("Fit Beta distribution to LAD...")
    # if (nrow(LAD_lim)>1){
    m <- fitdistrplus::fitdist(as.numeric(LAD$a)/90, 'beta', method='mle')
    
    # Get alpha and beta parametrs
    alpha0 <- m$estimate[1] # parameter 1
    beta0 <- m$estimate[2] # parameter 2
    beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01),
                                                     alpha0,beta0))
    param<-data.frame(alpha0,beta0)
    
    print("Finalizing TLSLeAF object...")
    
    TLSLeAF.dat<-new("TLSLeAF",
                     parameters=data.frame(c(file=input_file,
                                             center,
                                             scatterLim=85,
                                             SS=0.02,
                                             scales=c(0.1,0.5,0.75),
                                             voxRes=voxRes,
                                             superDF=TRUE)),
                     dat=as.data.frame(dat.class),
                     voxels=as.data.frame(voxels),
                     LAD=LAD,
                     Beta_parameters=param,
                     beta=beta,
                     G=G)
    
    # remove(dat)
    remove(dat.class)
    remove(voxels)
    remove(LAD)
    gc()
    return(TLSLeAF.dat)
  }
  
  if(clean){
    clean.temp(output_file,c2c.file, clean=TRUE)
  }
  
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

clean.temp<-function(output_file, c2c.file, clean=TRUE){
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


pgapF<-function (tls.scan.sub, row_res, col_res, z_res) {
  total<-list<-pgap<-NULL
  total<-(row_res/150)*col_res*5
  list <- seq(floor(min(tls.scan.sub[3])), 
              ceiling(max(tls.scan.sub[3])), 
              by = z_res)
  
  pgap <- NULL
  for (i in 1:length(list)) {
    pgap[i] <- 1-(length(subset(tls.scan.sub$z, tls.scan.sub$z<list[i])))/total
  }
  pgap<- (1-max(pgap))+pgap
  return(pgap)
}

a_vai<-function(tls.scan, a_ls, z_res, col_res, row_res){
  vai_a<-list()
  for (a in 1:length(a_ls)){
    
    tls.scan.sub<-subset(tls.scan,tls.scan$inc>a_ls[a]&tls.scan$inc<(a_ls[a]+5))
    list <- seq(floor(min(tls.scan.sub[3])), 
                ceiling(max(tls.scan.sub[3])), 
                by = z_res)
    
    pgap<-pgapF(tls.scan.sub, row_res, col_res, z_res)
    
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
      # tls.scan = gsub(".ptx","", files[l]),
      a = a_ls[a],
      z = f.prime.sub[,2],
      dx = f.prime.sub[,1])
    
    
  }
  return(data.frame(do.call(rbind,vai_a)))
}

vai<-function (tls.scan, a_range, a_bin, z_res, col_res, row_res, name){
  # tls.scan<-RZA(tls.scan)
  a_ls<-seq(min(a_range),max(a_range), by = a_bin) # zenith angle range
  
  vai_out<-data.frame(a_vai(tls.scan, a_ls, z_res, col_res, row_res),
                      tls.scan = name)
  return(vai_out)
}



vai_weighted <- function(a_vai, adj_z, z_res) {
  a_vai <- a_vai[a_vai$z>0,]
  a_vai$z_adj<-a_vai$z+adj_z
  a_vai$z_bin<-floor(a_vai$z_adj)
  #weight the profiles by zenith angle
  a_vai$weight<-sin(a_vai$a*(pi/180))
  a_vai_weight<-a_vai %>%
    group_by(tls.scan, z_bin) %>% 
    mutate(dxw = weighted.mean(dx, weight))
  return(a_vai_weight)
}


ptx.header<-function(input_file){
  list(
    name = input_file,
    col_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, header = FALSE)[1,]),
    row_res = as.numeric(read.csv(input_file, sep = "", skip = 0,nrows = 2, header = FALSE)[2,]),
    scan_center = as.matrix(read.csv(input_file, sep = "", skip = 2,nrows = 1, header = FALSE)),
    reg = as.matrix(read.csv(input_file, sep = "", skip = 3,nrows = 3, header = FALSE)),
    trans_matrix = as.matrix(read.csv(input_file, sep = "", skip = 6,nrows = 4, header = FALSE))
  )
}

pgap.angle<-function(tls.scan,
                     z_res=1,
                     adj_z=1.6,
                     a_range=c(10,60),
                     a_bin = 5,
                     plane = FALSE,
                     header,...){
  
  
  
  #calculate inclination angle
  inc<-RZA(tls.scan, deg=TRUE)
  tls.scan$inc<-inc
  remove(inc)
  gc()
  gc()
  
  # zenith angle range
  a_ls<-seq(min(a_range),max(a_range), by = a_bin) 
  
  pgap_ls<-list()
  
  pgap_ls<-lapply(a_ls, function(x){
    tls.scan.sub<-subset(tls.scan,tls.scan$inc>x&tls.scan$inc<(x+a_bin))
    tls.scan.sub<-na.exclude(tls.scan.sub)
    # if(nrow(tls.scan.sub) ==0) next
    list <- seq(floor(min(tls.scan.sub[3])), 
                ceiling(max(tls.scan.sub[3])), 
                by = z_res)
    
    pgap<-pgapF(tls.scan.sub, header$row_res, header$col_res, z_res)
    # plot(pgap)
    
    remove(tls.scan.sub)
    gc()
    
    return(data.frame(list,pgap, inc_bin=x))
  })
  remove(tls.scan)
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


G_calculations<-function(dat,
                         LAD,
                         voxels){
  
  G_ls<-list()
  if(file.exists(dat)){
    # dat<-NULL
    gc()
    dat<-fread(dat)
    dat<-dat[dat$class==0,]
    
    if(nrow(dat)>0){
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
      for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(dip_bin=G_test$dip_bin[G_test$dip_bin==c(0:90)[i]],
                                                         count=G_test$count[G_test$dip_bin==c(0:90)[i]],
                                                         p=G_test$p[G_test$dip_bin==c(0:90)[i]],
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
      
      # plot(G_calc$inc_bin,
      #      G_calc$G_p)
      
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
    } 
  }
  
  if(file.exists(voxels)){
    dat<-fread(voxels)
    if(nrow(dat)>0){
      
      dat$inc<-RZA(dat[,1:3], deg = TRUE)
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
      for(i in 1:length(0:90)) ran_G_ls[[i]]<-data.frame(dip_bin=G_test$dip_bin[G_test$dip_bin==c(0:90)[i]],
                                                         count=G_test$count[G_test$dip_bin==c(0:90)[i]],
                                                         p=G_test$p[G_test$dip_bin==c(0:90)[i]],
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
    }
  }
  
  fwrite(do.call(rbind, G_ls),
         file = gsub("angles_class_voxels.asc", "G_function.txt",voxels),
         sep= " ")
  dat<-NULL
  G_ls<-NULL
  gc()
  # return(do.call(rbind, G_ls))
}

