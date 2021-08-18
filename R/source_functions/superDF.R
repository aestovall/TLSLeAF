if (superDF){
    
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
  