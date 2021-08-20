TLSLeAF<-function(input_file,
                  overwrite=TRUE,
                  center, 
                  scatterLim=85,
                  SS=0.02, 
                  scales=c(0.1,0.5,0.75),
                  # rf_model,
                  rf_model_path,
                  correct.topography=TRUE,
                  voxRes = 0.1,
                  minVoxDensity=5,
                  superDF=TRUE,
                  clean=TRUE,...){
  
  # keep.vars<-ls()
  
  print(paste0("Processing ", input_file))
  
  #names of output files
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  gc()
  
  #Calculate normals for gridded TLS point cloud
  source("R/source_functions/normalCalc.R", local=TRUE)
  
  #Calculate scattering angle and leaf angle
  source("R/source_functions/angleCalcWrite.R", local=TRUE)
  
  #Classify wood and leaf from random forest classifier
  source("R/source_functions/classMetricCalc.R", local=TRUE)
  source("R/source_functions/rf_predict.R", local=TRUE)
  
  status<-TRUE
  print("Done")
  
  source("R/source_functions/voxel_beta_fit.R", local=TRUE)
  
  
  source("R/source_functions/G_calculations.R", local=TRUE)
  
  source("R/source_functions/superDF.R", local=TRUE)
  
  source("R/source_functions/clean.R", local=TRUE)
  
  if(superDF) return(TLSLeAF.dat)
  
  # plot(density(TLSLeAF.dat@LAD$a), main="LAD")
  # plot(TLSLeAF.dat@G[,1:2], main="G-functions")
  
}
