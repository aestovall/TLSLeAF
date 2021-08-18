#### EDIT THE SETUP FILE FOR YOUR SYSTEM ####

#Make sure the setup file is correct and save
# The most important step is to have Cloud Compare installed
# and add the executable into the setup file
file.edit('R/000_setup_windows.R')

# Run the setup file


files<-list.files("D:/PACE_TLS/Takeout/Drive/pace_TLS/input",pattern = "ptx", 
                  recursive = FALSE, 
                  full.names = TRUE)
source('R/000_setup_windows.R')

# remove(rf_model)


lapply(57:length(files), function(x){
  # source('R/000_setup_windows.R')
  
  rf_model<-readRDS("leaf_wood_class_RF.rds")
  
  ###input file ####
  
  input_file = files[x]
  print(input_file)
  
  #Find center coordinates. Modify according to file format
  center<-fread(input_file, nrows=4, header=FALSE)[1,]
  colnames(center)<-c("x","y","z")
  
  TLSLeAF.df<-NULL
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
                      clean=FALSE)
  TLSLeAF.df<-NULL
  remove(rf_model)
  
  # if(is.null(TLSLeAF.df)) next
  
  # header<- ptx.header(input_file)
  # 
  # #read in the scan file and rename columns
  # tls.scan<-fread(input_file, skip = 10, header = FALSE, select=c(1:3))
  # 
  # for(ii in 1:length(tls.scan)){
  #   if (class(tls.scan[[ii]])=="character"){
  #     tls.scan[[ii]]<- as.numeric( tls.scan[[ii]] )
  #   }
  # }
  # 
  # colnames(tls.scan)[1:3]<-c("x","y","z")
  # 
  # pgap_all<-pgap.angle(tls.scan, header=header)
  # remove(tls.scan)
  # gc()
  
  # PAVD<-NULL
  # G<-TLSLeAF.df@G
  
  # PAVD<-pgap2PAI(pgap_all = pgap_all,
  #                method = "C",
  #                G=G[G$assumption=="random"&
  #                      G$density=="normalized",])
  
  # write.csv(TLSLeAF.df@LAD, 
  #           gsub(".ptx","_LAD_voxels.csv", input_file), 
  #           row.names = FALSE)
  
  # write.csv(TLSLeAF.df@G, 
  #           gsub(".ptx","_G_function.csv", input_file), 
  #           row.names = FALSE)
  
  # write.csv(PAVD,
  #           gsub(".ptx","_PAVD.csv", input_file),
  #           row.names = FALSE)
  
  # output_file = gsub(".ptx",".asc", input_file)
  # angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  # c2c.file<-angle.file.name
  # class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  # 
  # clean.temp(output_file, c2c.file, clean=TRUE)
  
  # remove(TLSLeAF.df)
  gc()
  # .rs.restartR()
})


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




