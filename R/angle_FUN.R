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

LAvoxel<-function(dat, res = 0.1){
  las<-LAS(dat[,1:3])
  crs(las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  las@data$Z<-dat$z_cor
  las@data$dip_deg<-dat$dip_deg
  las@data$count<-1
  
  voxels<-grid_metrics3d(las, mean(dip_deg), res=res)
  colnames(voxels)[4]<-"dip_dir"
  voxels$dip_dir_sd<-grid_metrics3d(las, sd(dip_deg), res=res)[,4]
  voxels$n<-grid_metrics3d(las, sum(count), res=res)[,4]
  
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

