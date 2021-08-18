
# file.edit('R/000_setup_windows.R')

files<-list.files("D:/PACE_TLS/Takeout/Drive/pace_TLS/input",pattern = "ptx", 
                  recursive = FALSE, 
                  full.names = TRUE)

# Run the setup file
source('R/000_setup_windows.R')


input_file
overwrite=FALSE
scatterLim=80
SS=0.02
scales=c(0.1,0.5,0.75)
rf_model=rf_model
# rf_model_path=rf_model_path,
correct.topography = TRUE
voxRes=0.1
minVoxDensity=5
superDF=TRUE
clean=TRUE

setDTthreads(threads = 1)
getDTthreads(verbose=TRUE)



for(i in 49:length(files)){
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
  
  print(paste0("Processing ", input_file))
  
  #names of output files
  output_file = gsub(".ptx",".asc", input_file)
  angle.file.name<-gsub(".ptx","_angles.asc", input_file)
  c2c.file<-angle.file.name
  class.file.name<-gsub(".ptx","_angles_class.asc", input_file)
  gc()
  
  print("Calculate normals for gridded TLS point cloud...")
  #Calculate normals for gridded TLS point cloud
  if(!file.exists(output_file)|
     overwrite) normalCalc(input_file)
  print("Done")
  
  while(!file.exists(output_file)) Sys.sleep(10)
  
  print("Calculate scattering angle and leaf angle...")
  #Calculate scattering angle and leaf angle
  if((length(rstudioapi::terminalBusy(rstudioapi::terminalList())[
    rstudioapi::terminalBusy(rstudioapi::terminalList())])==0)&
    (!file.exists(angle.file.name)|
     overwrite)){
    # remove(dat)
    # gc()
    
    dat<-readTLSnorms(output_file)
    
    dat.angle<-angleCalc(dat,
                         center,
                         scatterLim)
    
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
  
  if(file.exists(c2c.file)&
     !file.exists(gsub(".asc","_rf_prep.asc",c2c.file))) {
    
    print("Prepping for leaf and wood classification...")
    rfPrep(c2c.file)
    
    print("Done")  
  }
  
  while(!file.exists(gsub(".asc","_rf_prep.asc",c2c.file))) Sys.sleep(10)
  
  if(!file.exists(class.file.name)&
     file.exists(gsub(".asc","_rf_prep.asc",c2c.file))) {
    
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
  
  # # while(!file.exists(gsub(".asc", paste("_Beta_distribution_alpha_", 
  # #                                       round(param[1],2),"_beta_",round(param[2],2),
  # #                                       ".txt", sep=""),class.file.name))) Sys.sleep(5)
  # 
  # if (superDF){
  #   
  #   dat.class<-LAD<-voxels<-G<-m<-beta<-NULL
  #   
  #   dat.class<-fread(class.file.name)
  #   LAD<-fread(gsub(".asc", "_LAD.txt",class.file.name))
  #   voxels<-fread(gsub(".asc", "_voxels.asc",class.file.name))
  #   G<-fread(gsub("angles_class_voxels.asc", "G_function.txt",gsub(".asc", "_voxels.asc",class.file.name)))
  #   
  #   # get_beta<-function(x){
  #   m<-NULL
  #   print("Fit Beta distribution to LAD...")
  #   # if (nrow(LAD_lim)>1){
  #   m <- fitdistrplus::fitdist(as.numeric(LAD$a)/90, 'beta', method='mle')
  #   
  #   # Get alpha and beta parametrs
  #   alpha0 <- m$estimate[1] # parameter 1
  #   beta0 <- m$estimate[2] # parameter 2
  #   beta<-data.frame(a= seq(0.01,0.98, 0.01),y=dbeta(seq(0.01,0.98, 0.01),
  #                                                    alpha0,beta0))
  #   param<-data.frame(alpha0,beta0)
  #   
  #   print("Finalizing TLSLeAF object...")
  #   
  #   TLSLeAF.dat<-new("TLSLeAF",
  #                    parameters=data.frame(c(file=input_file,
  #                                            center,
  #                                            scatterLim=85,
  #                                            SS=0.02,
  #                                            scales=c(0.1,0.5,0.75),
  #                                            voxRes=voxRes,
  #                                            superDF=TRUE)),
  #                    dat=as.data.frame(dat.class),
  #                    voxels=as.data.frame(voxels),
  #                    LAD=LAD,
  #                    Beta_parameters=param,
  #                    beta=beta,
  #                    G=G)
  #   
  #   # remove(dat)
  #   remove(dat.class)
  #   remove(voxels)
  #   remove(LAD)
  #   gc()
  #   return(TLSLeAF.dat)
  
  if(clean){
    clean.temp(output_file,c2c.file, clean=TRUE)
  }
  
}

