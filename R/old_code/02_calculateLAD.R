# Calculate LAD

#Create a surface model and topographically correct the TLS data
if(correct.topography == TRUE){
  las<-LAS(dat[,1:3])
  crs(las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  
  r <- raster(xmn=-200, xmx=200, ymn=-200, ymx=200, resolution = 5)
  topo<-grid_metrics(las, quantile(Z, 0.01), r)
  plot(topo, col = viridis::viridis(250))
  
  crs(topo)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  slope<-terrain(topo, opt = "slope", unit = "degrees", neighbors = 8)
  plot(slope)
  
  topo[slope>40]<-NA
  setMinMax(topo)
  
  topo.df<-as.data.frame(rasterToPoints(topo))
  colnames(topo.df)<-c("X","Y","Z")
  
  ws <- seq(3,12, 3)
  th <- seq(0.1, 1.5, length.out = length(ws))
  
  topo.las<-LAS(topo.df)
  crs(topo.las)<-"+proj=eqc +lat_ts=0 +lat_0=0 +lon_0=0 +x_0=0 +y_0=0 +a=6371007 +b=6371007 +units=m +no_defs"
  
  ground<-lasground(topo.las, pmf(ws, th), last_returns = FALSE)
  
  topo.las.r<-grid_terrain(ground, res = 5, knnidw(k = 21))
  plot(topo.las.r)
  
  las<- las - topo.las.r
  
  dat$z_cor<-las@data$Z
  dat[dat$z_cor<1|dat$z_cor>70,]<-NA
} else dat$z_cor<-dat$Z

#Leaf angle distribution
LAD<-hist(dat$dip_deg, breaks = c(0:91))
LAD<-data.frame(ID = gsub(".asc","",input_file),
                a = LAD$breaks[-91], 
                count = LAD$counts,
                density = LAD$density)

plot(LAD$a,LAD$density, col = "white")
lines(LAD$a,LAD$density, col = "blue")

  