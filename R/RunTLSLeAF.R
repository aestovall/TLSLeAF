#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file
file.edit('R/000_setup_windows.R')

# Run the setup file


files<-list.files("../TLS_Leaf_Angle/input",pattern = "ptx", 
                  recursive = FALSE, 
                  full.names = TRUE)
source('R/000_setup_windows.R')

for(i in 46:length(files)){
  # source('R/000_setup_windows.R')
  
  gc()
  gc()
  
  ###input file ####
  
  input_file = files[i]
  
  #Find center coordinates. Modify according to file format
  center<-fread(input_file, nrows=4, header=FALSE)[1,]
  colnames(center)<-c("x","y","z")
  
  # TLSLeAF.df<-NULL
  TLSLeAF.df<-TLSLeAF(input_file, 
                      overwrite=FALSE,
                      center, 
                      scatterLim=80,
                      SS=0.02, 
                      scales=c(0.1,0.5,0.75),
                      rf_model,
                      # rf_model_path=rf_model_path,
                      correct.topography = TRUE,
                      voxRes=0.1,
                      minVoxDensity=5,
                      superDF=TRUE,
                      clean=TRUE)
  
  # if(is.null(TLSLeAF.df)) next
  
  header<- ptx.header(input_file)

  #read in the scan file and rename columns
  tls.scan<-fread(input_file, skip = 10, header = FALSE, select=c(1:3))

  for(ii in 1:length(tls.scan)){
    if (class(tls.scan[[ii]])=="character"){
      tls.scan[[ii]]<- as.numeric( tls.scan[[ii]] )
    }
  }

  colnames(tls.scan)[1:3]<-c("x","y","z")

  pgap_all<-pgap.angle(tls.scan, header=header)
  remove(tls.scan)
  gc()

  PAVD<-NULL
  G<-TLSLeAF.df@G
  
  PAVD<-pgap2PAI(pgap_all = pgap_all,
                 method = "C",
                 G=G[G$assumption=="random"&
                       G$density=="normalized",])
  # 
  # PAVD.sum<-aggregate(f.pred.prime~list,FUN="mean",PAVD)
  # PAVD.sum$sd<-aggregate(f.pred.prime~list,FUN="sd",PAVD)[,2]
  # 
  # PAVD.sum<-aggregate(f.prime~list,FUN="mean",PAVD)
  # PAVD.sum$sd<-aggregate(f.prime~list,FUN="sd",PAVD)[,2]
  # 
  # # library(ggplot2)
  # # ggplot(PAVD.sum, aes(x= list, y= f.pred.prime))+
  # #   geom_ribbon(aes(ymin=f.pred.prime-sd, ymax=f.pred.prime+sd), 
  # #               alpha=0.5, fill="forestgreen", color=NA)+
  # #   geom_path()+
  # #   coord_flip()+
  # #   xlab("Height (m)")+
  # #   ylab(bquote(PAVD~(m^2~m^{-2})))
  # # ggplot(PAVD.sum, aes(x= list, y= f.prime))+
  # #   geom_ribbon(aes(ymin=f.prime-sd, ymax=f.prime+sd), 
  # #               alpha=0.5, fill="forestgreen", color=NA)+
  # #   geom_path()+
  # #   coord_flip()+
  # #   xlab("Height (m)")+
  # #   ylab(bquote(PAVD~(m^2~m^{-2})))
  # 
  # plot(PAVD.sum$list, PAVD.sum$f.pred.prime)
  # lines(PAVD.sum$list, PAVD.sum$f.pred.prime)
  # 
  write.csv(TLSLeAF.df@LAD, 
            gsub(".ptx","_LAD_voxels.csv", input_file), 
            row.names = FALSE)
  
  write.csv(TLSLeAF.df@G, 
            gsub(".ptx","_G_function.csv", input_file), 
            row.names = FALSE)
  
  write.csv(PAVD,
            gsub(".ptx","_PAVD.csv", input_file),
            row.names = FALSE)
  
  # par(mfrow=c(2,2))
  # plot(df@beta)
  # hist(df@LAD$a.a, breaks=90)
  # plot(aggregate(a.a~a.Z, FUN="mean", df@LAD))
  # plot(G$inc_bin,G$G_p, col="black")
  
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  
  clean.temp(output_file, c2c.file, clean=TRUE)
  
  # env.vars<-ls()
  # remove(list=env.vars[env.vars!="files"])
  remove(TLSLeAF.df)
  gc()
  
  # Sys.sleep(20)
  
}
  

library(parallel)

cl<-makeCluster(getOption("cl.cores", 4))


parLapply(cl, files[1:4], function(x){
  source('R/000_setup_windows.R')
  
  ###input file ####
  
  # input_file = files[i]
  input_file = x
  
  #Find center coordinates. Modify according to file format
  center<-fread(input_file, nrows=4, header=FALSE)[1,]
  colnames(center)<-c("x","y","z")
  
  df<-TLSLeAF(input_file, 
              overwrite=FALSE,
              center, 
              scatterLim=80,
              SS=0.02, 
              scales=c(0.1,0.5,0.75),
              rf_model,
              correct.topography = TRUE,
              voxRes=0.1,
              minVoxDensity=5,
              superDF=TRUE)
  
  # plot(df@beta)
  # hist(df@LAD$a.a, breaks=90)
  # plot(aggregate(a.a~a.Z, FUN="mean", df@LAD))
  
  G<-G_calculations(df)
  
  write.csv(gsub(".ptx","_G_function.csv", input_file), 
            row.names = FALSE)
  
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  
  clean.temp(output_file, c2c.file)
  
  gc()
  gc()
  
})
stopCluster(cl)




