file.edit('R/000_setup_windows.R')

# Run the setup file


files<-list.files("D:/PACE_TLS/Takeout/Drive/pace_TLS/input",pattern = "ptx", 
                  recursive = FALSE, 
                  full.names = TRUE)
source('R/000_setup_windows.R')

for(i in 2:length(files)){
  # source('R/000_setup_windows.R')
  print(i)
  
  rf_model<-readRDS("leaf_wood_class_RF.rds")
  
  gc()
  gc()
  
  ###input file ####
  input_file = files[i]
  print(input_file)
  #Find center coordinates. Modify according to file format
  center<-fread(input_file, nrows=4, header=FALSE)[1,]
  colnames(center)<-c("x","y","z")
  
  TLSLeAF(input_file, 
          overwrite=TRUE,
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
  
}
