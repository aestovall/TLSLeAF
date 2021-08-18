file.edit('R/000_setup_windows.R')

files<-list.files("D:/PACE_TLS/Takeout/Drive/pace_TLS/input",pattern = "ptx", 
                  recursive = FALSE, 
                  full.names = TRUE)

# Run the setup file
source('R/000_setup_windows.R')

for(i in 42:length(files)){
  # source('R/000_setup_windows.R')
  print(i)
  
  rf_model<-readRDS("leaf_wood_class_RF.rds")
  
  ###input file ####
  input_file = files[i]
  print(input_file)
  
  #Find center coordinates. Modify according to file format
  center<-fread(input_file, nrows=4, header=FALSE)[1,]
  colnames(center)<-c("x","y","z")
  
  TLSLeAF.df<-NULL
  gc()
  TLSLeAF.df<-TLSLeAF(input_file,
                      overwrite=FALSE,
                      center,
                      scatterLim=80,
                      SS=0.02,
                      scales=c(0.1,0.5,0.75),
                      rf_model=rf_model,
                      # rf_model_path=rf_model_path,
                      correct.topography = TRUE,
                      voxRes=0.1,
                      minVoxDensity=5,
                      superDF=TRUE,
                      clean=TRUE)
  
  gc()
  
  # output_file = gsub(".ptx",".asc", input_file)
  # angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  # c2c.file<-angle.file.name
  # class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  # 
  # clean.temp(output_file, c2c.file, clean=TRUE)
  
  gc()
}
